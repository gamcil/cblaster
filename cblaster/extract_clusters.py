#!/usr/bin/env python3

"""Extract clusters from session files into genbank or fasta files"""


import logging
import os
import time
import requests
import json
from g2j.classes import Organism


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation


from cblaster.classes import Session
from cblaster.extract import organism_matches, parse_organisms, parse_scaffolds
from cblaster.helpers import efetch_sequences


LOG = logging.getLogger(__name__)

# from https://www.ncbi.nlm.nih.gov/books/NBK25497/
MAX_REQUEST_SIZE = 500
MIN_TIME_BETWEEN_REQUEST = 0.34  # seconds


def parse_numbers(cluster_numbers):
    """parse clutser numbers from user input

    Args:
        cluster_numbers (list): a list of numbers or ranges of numbers
    Returns:
        list of integer numbers
    """
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
        Set of tuples in the form (cblaster.Cluster object, scaffold_accession of cluster, organism_name of cluster)
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
    """Get all clusters of an organism

    Args:
        organism (cblaster.Organism): cblaster Organism object
    Returns:
        list of tuples in form (cblaster.Cluster object, scaffold_accession of cluster, organism_name of cluster)
    """
    selected_clusters = []
    for scaffold in organism.scaffolds.values():
        selected_clusters.extend(get_scaffold_clusters(scaffold, organism.name))
    return selected_clusters


def get_scaffold_clusters(scaffold, organism_name, start=None, end=None):
    """Get all clusters on a scaffold

    Args:
        scaffold (cblaster.Scaffold): cblaster scaffold object
        organism_name(str): name of the organism the scaffold is from
        start (int): start coordinate of cluster must be bigger then start
        end (int): end coordiante of cluster must be smaller then stop
    Returns:
        list of tuples in form (cblaster.Cluster object, scaffold_accession of cluster, organism_name of cluster)
    """
    selected_clusters = []
    for cluster in scaffold.clusters:
        # check if the cluster is within the range given for the scaffold
        if (start and start >= cluster.start) or (end and end <= cluster.end):
            continue
        selected_clusters.append((cluster, scaffold.accession, organism_name))
    return selected_clusters


def create_genbanks_from_clusters(session, cluster_hierarchy, output_dir, prefix):
    """Create genbank files for each selected cluster

    Args:
        session (Session): a cblaster session object
        cluster_hierarchy (Set): a set of selected clusters
        output_dir (string): path to a directory for writing the output files
        prefix (string): string to start the file name of each cluster with
    """
    if session.params["mode"] == "remote":
        protein_sequences = efetch_protein_sequences(cluster_hierarchy)
        nucleotide_sequences = efetch_nucleotide_sequence(cluster_hierarchy)
    elif session.params["mode"] == "local":
        protein_sequences, nucleotide_sequences = database_fetch_sequences(session.params["json_db"], cluster_hierarchy)
    else:
        raise NotImplementedError(f"No protocol for mode {session.params['mode']}")

    # generate genbank files for all the required clusters
    for cluster, scaffold_accession, organism_name in cluster_hierarchy:
        cluster_prot_sequences = {subject.name: protein_sequences[subject.name] for subject in cluster.subjects}
        cluster_nuc_sequence = nucleotide_sequences[scaffold_accession]
        with open(f"{output_dir}{os.sep}{prefix}cluster{cluster.number}.gb", "w") as f:
            record = cluster_to_record(cluster, cluster_prot_sequences, cluster_nuc_sequence, organism_name, scaffold_accession)
            SeqIO.write(record, f, 'genbank')
        LOG.debug(f"Created {prefix}cluster{cluster.number}.gb file for cluster {cluster.number}")


def database_fetch_sequences(json_db, cluster_hierarchy):
    """Fetch sequences from an offline json_db

    Args:
        json_db (str): path to a json_db file
        cluster_hierarchy (Set): set of tuples in the form (cblaster.Cluster object,
         scaffold_accession of cluster, organism_name of cluster)
    Returns:
        tuple of a dict of protein sequences keyed on protein_names and a dict of
         nucleotide sequences keyed on scaffold_accessions
    """
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
                            break
    return prot_sequences, nuc_sequences


def efetch_protein_sequences(cluster_hierarchy):
    """eFetch all protein sequences for all the clusters in the cluster_hierarchy

    Args:
        cluster_hierarchy (Set): set of tuples in the form (cblaster.Cluster object,
         scaffold_accession of cluster, organism_name of cluster)
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
    """eFetch all nucleotide sequences for all the clusters in the cluster_hierarchy

    Args:
        cluster_hierarchy (Set): set of tuples in the form (cblaster.Cluster object,
         scaffold_accession of cluster, organism_name of cluster)
    Returns:
        a dictionary with nucleotide sequences keyed on scaffold accession
    """
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
        LOG.info(f"Querying NCBI for sequence for {scaffold_accession} from {cluster.start} to {cluster.end}")
        LOG.debug(f"Efetch URL: {response.url}")

        if response.status_code != 200:
            raise requests.HTTPError(
                f"Error fetching sequences from NCBI [code {response.status_code}]."
                " Incorect scaffold accession?"
            )
        # only save the sequence not the fasta header
        sequences[scaffold_accession] = response.text.split("\n", 1)[1].replace("\n", "")

        passed_time = time.time() - start_time
    return sequences


def cluster_to_record(cluster, cluster_prot_sequences, cluster_nuc_sequence, organism_name, scaffold_accession):
    """Convert a cblaster.Cluster object into a Bio.Seqrecord object

    Args:
        cluster (cblaster.Cluster): cblasdter Cluster object
        cluster_prot_sequences (dict): dictionary of protein sequences keyed on protein ids
        cluster_nuc_sequence (dict): dictionary of nucleotide sequences keyed on scaffold accession
        organism_name (str): name of the organism the cluster is originated from
        scaffold_accession (str): accession of the scaffold the cluster belongs to
    Returns:
        a Bio.Seqrecord object
    """
    nuc_seq_obj = Seq(cluster_nuc_sequence, IUPAC.unambiguous_dna)

    # create the record
    record = SeqRecord(nuc_seq_obj,
                       id=scaffold_accession,
                       name=scaffold_accession,
                       description=f"Genes for cluster {cluster.number} on scaffold {scaffold_accession}")
    feature = SeqFeature(FeatureLocation(start=cluster.start, end=cluster.end), type="SOURCE",
                         qualifiers={"organism": organism_name})
    record.features.append(feature)

    subjects = {subject.name: subject for subject in cluster.subjects}
    for name, sequence in cluster_prot_sequences.items():
        subject = subjects[name]
        cds_feature = SeqFeature(FeatureLocation(start=subject.start - cluster.start, end=subject.end- cluster.start),
                                 type="CDS", qualifiers={"protein_id": subject.name, "translation": sequence})
        record.features.append(cds_feature)
    return record


def extract_clusters(
    session,
    output_dir,
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

    LOG.info(f"Writing genbank files")
    create_genbanks_from_clusters(session, cluster_hierarchy, output_dir, prefix)

    LOG.info(f"All clusters have been written to files. Output can be found at {output_dir}")
    LOG.info("Done!")