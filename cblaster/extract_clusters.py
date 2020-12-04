#!/usr/bin/env python3

"""Extract clusters from session files into genbank or fasta files"""


import logging
import os
import time
import requests
import json
from g2j.classes import Organism


from cblaster.classes import Session
from cblaster.extract import organism_matches, parse_organisms, parse_scaffolds
from cblaster.helpers import efetch_sequences


LOG = logging.getLogger(__name__)

# from https://www.ncbi.nlm.nih.gov/books/NBK25497/
MAX_REQUEST_SIZE = 500
MIN_TIME_BETWEEN_REQUEST = 0.34  # seconds

# constants for genbank files
MAX_LINE_LENGTH = 80
FEATURE_SPACE_LENGTH = 21
NUCLEOTIDE_LINE_LENGTH = 60


def parse_numbers(cluster_numbers):
    chosen_cluster_numbers = set()
    for number in cluster_numbers:
        try:
            if "-" in number:
                start, stop = number.split("-")
                chosen_cluster_numbers.update(range(int(start), int(stop) + 1))
            else:
                chosen_cluster_numbers.add(int(number))
        except ValueError:
            LOG.warning(f"Cannot extract cluster '{number}' the given number is not a valid integer")
    return chosen_cluster_numbers


def extract_cluster_hierarchies(
    session,
    cluster_numbers,
    score_threshold,
    organisms,
    scaffolds
):
    """Filter out selected clusters
    Args:
        session (Session): A session object
        cluster_numbers (list): Numbers of clusters or a range of numbers eg 1-5
        score_threshold (float): Minimum score a cluster needs to have in order to be included
        organisms (list): Regex patterns for organisms of which all clusters need to be extracted
        scaffolds (list): Names of scaffolds of which all clusters need to be extracted
    Returns:
        Set of cluster names e.g. posive integers
    """
    # no filter options return all clusters
    selected_clusters = set()
    if not cluster_numbers and not score_threshold and not organisms and not scaffolds:
        for organism in session.organisms:
            selected_clusters.update(get_organism_clusters(organism))
        return selected_clusters

    # prepare the filters defined by the user

    if cluster_numbers:
        cluster_numbers = parse_numbers(cluster_numbers)
    if organisms:
        organisms = parse_organisms(organisms)
    if scaffolds:
        scaffolds = parse_scaffolds(scaffolds)

    # actually filter out the clusters
    for organism in session.organisms:
        if organisms and organism_matches(organism.name, organisms):
            selected_clusters.update(get_organism_clusters(organism))
            continue
        for scaffold in organism.scaffolds.values():
            if scaffolds and scaffold.accession in scaffolds:
                selected_clusters.update(get_scaffold_clusters(scaffold, organism.name, **scaffolds[scaffold.accession]))
                continue
            for cluster in scaffold.clusters:
                if cluster_numbers and cluster.number in cluster_numbers or \
                        (score_threshold and cluster.score >= score_threshold):
                    selected_clusters.add((cluster, scaffold.accession, organism.name))
    return selected_clusters


def get_organism_clusters(organism):
    """Select all clusters of an organism"""
    selected_clusters = []
    for scaffold in organism.scaffolds.values():
        selected_clusters.extend(get_scaffold_clusters(scaffold, organism.name))
    return selected_clusters


def get_scaffold_clusters(scaffold, organism_name, start=None, end=None):
    """Select all clusters on a scaffold"""
    selected_clusters = []
    for cluster in scaffold.clusters:
        # check if the cluster is within the range given for the scaffold
        if (start and start >= cluster.start) or (end and end <= cluster.end):
            continue
        selected_clusters.append((cluster, scaffold.accession, organism_name))
    return selected_clusters


def create_files_from_clusters(session, cluster_hierarchy, output_dir, prefix, file_format):
    """Create file_fromat files for each selected cluster
    Args:
        session (Session): a cblaster session object
        cluster_hierarchy (Set): a set of selected clusters
        output_dir (string): path to a directory for writing the output files
        prefix (string): string to start the file name of each cluster with
        file_format (string):  the format of the output cluster files either
         genbank or fasta
    """
    if session.params["mode"] == "remote":
        protein_sequences = efetch_protein_sequences(session)
        nucleotide_sequences = efetch_nucleotide_sequence(cluster_hierarchy)
    elif session.params["mode"] == "local":
        protein_sequences, nucleotide_sequences = database_fetch_sequences(session.params["json_db"], cluster_hierarchy)
    else:
        raise NotImplementedError(f"No protocol for mode {session.params['mode']}")

    # generate genbank files for all the required clusters
    for cluster, scaffold_accession, organism_name in cluster_hierarchy:
        cluster_prot_sequences = {subject.name: protein_sequences[subject.name] for subject in cluster.subjects}
        cluster_nuc_sequence = nucleotide_sequences[scaffold_accession]
        if file_format == "genbank":
            with open(f"{output_dir}{os.sep}{prefix}cluster{cluster.number}.gb", "w") as f:
                f.write(cluster_to_genbank(cluster, cluster_prot_sequences, cluster_nuc_sequence, organism_name, scaffold_accession))
        elif file_format == "fasta":
            with open(f"{output_dir}{os.sep}{prefix}cluster{cluster.number}.fa", "w") as f:
                f.write(cluster_to_fasta(cluster, cluster_prot_sequences, organism_name, scaffold_accession))
        else:
            # this should not be possible
            raise NotImplementedError(f"File format {file_format} is not supported.")
        LOG.debug(f"Created {file_format} file for cluster {cluster.number}")


def database_fetch_sequences(json_db, cluster_hierarchy):

    # read the json database
    with open(json_db, "r") as f:
        # g2j Organism objects not cblaster Organism objects
        organism_objs = [Organism.from_dict(dct) for dct in json.load(f)]

    # extract all sequences from the database
    organisms_dict = {organism.name: organism for organism in organism_objs}
    identifiers = ("protein_id", "locus_tag", "gene", "ID", "Name", "label")
    prot_sequences = dict()
    nuc_sequences = dict()
    for cluster, scaffold_accession, organism_name in cluster_hierarchy:
        organism = organisms_dict[organism_name]
        for scaffold in organism.scaffolds:
            if scaffold.accession != scaffold_accession:
                continue
            nuc_sequences[scaffold_accession] = scaffold.sequence[cluster.start : cluster.end + 1]
            for feature in scaffold.features:
                if feature.type != "CDS":
                    continue
                for subject in cluster.subjects:
                    for identifier in identifiers:
                        if identifier in feature.qualifiers and subject.name == feature.qualifiers[identifier]:
                            prot_sequences[subject.name] = feature.qualifiers["translation"]
    return prot_sequences, nuc_sequences


def efetch_protein_sequences(cluster_hierarchy):
    """Fetch all protein sequences for all the clusters in cluster_numbers
    Args:
        session (Session): a cblaster session object
        cluster_hierarchy (Set): a set of selected clusters
    Returns:
        a dictionary with protein sequences keyed on protein names
    """
    # first collect all names to do the fetching all at once
    sequence_names = []
    for cluster, sc_name, or_name in cluster_hierarchy:
        sequence_names.extend([subject.name for subject in cluster.subjects])

    # request sequences in batces of MAX_REQUEST_SIZE every 0.34 seconds (no more then 3 requests per second)
    sequences = dict()
    passed_time = 0
    for i in range(int(len(sequence_names) / MAX_REQUEST_SIZE) + 1):
        if passed_time < MIN_TIME_BETWEEN_REQUEST:
            time.sleep(MIN_TIME_BETWEEN_REQUEST - passed_time)
        start_time = time.time()
        subset_sequence_names = sequence_names[i * MAX_REQUEST_SIZE: (i + 1) * MAX_REQUEST_SIZE]
        sequences.update(efetch_sequences(subset_sequence_names))
        passed_time = time.time() - start_time
    return sequences


def efetch_nucleotide_sequence(cluster_hierarchy):
    sequences = dict()
    passed_time = 0
    for cluster, scaffold_accession, org_name in cluster_hierarchy:
        if passed_time < MIN_TIME_BETWEEN_REQUEST:
            time.sleep(MIN_TIME_BETWEEN_REQUEST - passed_time)
        start_time = time.time()
        response = requests.post(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
            params={"db": "sequences", "rettype": "fasta",
                    "seq_start": str(cluster.start), "seq_stop": str(cluster.end),
                    "strand": "1"},
            files={"id": scaffold_accession}
        )
        LOG.debug(f"Efetch sequence for {scaffold_accession} from {cluster.start} to {cluster.end}")
        LOG.debug(f"Efetch URL: {response.url}")

        if response.status_code != 200:
            raise requests.HTTPError(
                f"Error fetching sequences from NCBI [code {response.status_code}]."
                " Incorect scaffold accession?"
            )
        # only save the sequence not the fasta header
        sequences[scaffold_accession] = response.text.split("\n", 1)[1]

        passed_time = time.time() - start_time
    return sequences


def cluster_to_genbank(cluster, cluster_prot_sequences, cluster_nuc_sequence, organism_name, scaffold_accession):
    """Write a Cluster object into a genbank string
    Args:
        cluster (Cluster): a cblaster Cluster object
        cluster_prot_sequences (dict): a dictionary linking sequence names to the aa
         sequence
        organism_name (string): name of the organism the cluster originated from
        scaffold_accession (string): accession number of the scaffold the organism
         originates from
    Returns:
        a genbank string
    """

    # first construct a simple header
    genbank_string = f"LOCUS       {scaffold_accession}               {cluster.end - cluster.start} bp    DNA\n"
    genbank_string += f"DEFINITION  Genes for cluster {cluster.number} on scaffold {scaffold_accession}\n"
    genbank_string += f"ACCESSION   {scaffold_accession}\n"
    genbank_string += f"VERSION     {scaffold_accession}.1\n"
    genbank_string += "KEYWORDS    \n"
    genbank_string += f"SOURCE      {organism_name}\n"
    genbank_string += f"  ORGANISM  {organism_name}\n"
    genbank_string += f"REFERENCE   1  (bases {cluster.start} to {cluster.end})\n"
    genbank_string += f"  AUTHORS   Creators of cblaster\n"
    # TODO: add a proper title
    genbank_string += f"  TITLE     \n"
    genbank_string += "FEATURES             Location/Qualifiers\n"
    genbank_string += f"     source          {cluster.start}..{cluster.end}\n"
    genbank_string += f"                     /organism=\"{organism_name}\"\n"

    subjects = {subject.name: subject for subject in cluster.subjects}
    # construct a CDS feature for every subject in the cluster
    for name, sequence in cluster_prot_sequences.items():
        subject = subjects[name]
        hit_name = max(subject.hits, key=lambda x: x.bitscore).subject
        if subject.strand == "-":
            location_string = f"complement({subject.start - cluster.start}..{subject.end - cluster.start})"
        else:
            location_string = f"{subject.start - cluster.start}..{subject.end - cluster.start}"
        genbank_string += f"     CDS             {location_string}\n"
        genbank_string += f"                     /protein_id=\"{hit_name}\"\n"
        genbank_string += collapse_translation(f"                     /translation=\"{sequence}\"")
    genbank_string += "ORIGIN      \n"

    genbank_string += collapse_nucleotide_sequence(cluster_nuc_sequence)
    genbank_string += "//\n"
    return genbank_string


def collapse_translation(line):
    """Format a translation string to be in genbank form
    e.g.
        /translation="MRTTAQATSWLRRYRPRPAATWRLVCFPYA
        VCFPYASGNATFYRQWAVRLPAEVEVVAVQYPGRLDRIHEPCVR
        ASGNATF"
    Args:
        line (string): aa string to format
    Returns:
        the same string formatted to be in genbank form
    """
    collapsed_translation = line[0:MAX_LINE_LENGTH] + "\n"
    aas_per_line = MAX_LINE_LENGTH - FEATURE_SPACE_LENGTH
    for i in range(aas_per_line, len(line), aas_per_line):
        collapsed_translation += " " * FEATURE_SPACE_LENGTH + line[i:i + aas_per_line] + "\n"
    return collapsed_translation


def collapse_nucleotide_sequence(line):
    """Format a nuccleotide sequence string to be in genbank form
    e.g.
         1 aatgggcaaa aagc... etc
        61 aaggtacccg ttta... etc
    Args:
        line (string): nucleotide string to format
    Returns:
        the same string formatted to be in genbank form
    """
    collapesed_sequence = ""
    line = line.replace("\n", "")
    for i in range(0, len(line), NUCLEOTIDE_LINE_LENGTH):

        sequence_line = line[i:i + NUCLEOTIDE_LINE_LENGTH].lower()
        # add a space every 10 characters to resemble the NCBI genbank format
        spaced_sequence_line = ' '.join(sequence_line[i:i + 10] for i in range(0, len(sequence_line), 10))
        collapesed_sequence += f"{i + 1:>10} {spaced_sequence_line}\n"
    return collapesed_sequence


def cluster_to_fasta(cluster, cluster_sequences, organism_name, scaffold_accession):
    """Make a fasta string from a cluster object
    Args:
        cluster (Cluster): a cblaster Cluster object
        cluster_sequences (dict): a dictionary linking sequence names to the aa
         sequence
        organism_name (string): name of the organism the cluster originated from
        scaffold_accession (string): accession number of the scaffold the organism
         originates from
    Returns:
        a string in fasta format
    """
    fasta_string = ""
    for name, sequence in cluster_sequences.items():
        fasta_string += f">{name} on Cluster={cluster.number}, Organism={organism_name}," \
            f" Scaffold={scaffold_accession}\n"
        fasta_string += collapse_protein_sequence(sequence)
    return fasta_string


def collapse_protein_sequence(line):
    """Format a line of any lenght to MAX_LINE_LENGHT
    Args:
        line (string): an amino acid protein sequence
    Returns:
        the same string with a newline every MAX_LINE_LENGHT amount of characters
    """
    collapsed_sequence = ""
    for i in range(0, len(line), MAX_LINE_LENGTH):
        collapsed_sequence += line[i:i + MAX_LINE_LENGTH] + "\n"
    return collapsed_sequence


def extract_clusters(
    session,
    output_dir,
    file_format="genbank",
    prefix="",
    cluster_numbers=None,
    score_threshold=None,
    organisms=None,
    scaffolds=None,
):
    """Extract Cluster objects from a Session object and write them to a file in a
    specified format
    Args:
        session (string): path to a session.json file
        output_dir (string): path to a directory for writing the output files
        file_format (string):  the format of the output cluster files either
         genbank or fasta
        prefix (string): string to start the file name of each cluster with
        cluster_numbers (list): cluster numbers to include
        score_threshold (float): minum score in order for a cluster to be included
        organisms (list): Organism filtering regular expressions, clusters for
         these organisms are included
        scaffolds(list): clusters on these scaffolds are included
    """
    LOG.info("Starting cblaster plotting of clusters using clinker")
    LOG.info("Loading session from: %s", session)
    with open(session) as fp:
        session = Session.from_json(fp)

    LOG.info("Extracting clusters that match the filters")
    cluster_hierarchy = extract_cluster_hierarchies(session, cluster_numbers, score_threshold, organisms, scaffolds)
    LOG.info(f"Extracted {len(cluster_hierarchy)} clusters.")
    if len(cluster_hierarchy) == 0:
        LOG.info("There are no clusters that meet the filtering criteria. Exiting...")
        raise SystemExit

    LOG.info(f"Writing {file_format} files")
    create_files_from_clusters(session, cluster_hierarchy, output_dir, prefix, file_format)

    LOG.info(f"All clusters have been written to files. Output can be found at {output_dir}")
    LOG.info("Done!")

