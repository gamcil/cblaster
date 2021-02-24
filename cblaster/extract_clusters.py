#!/usr/bin/env python3

"""Extract clusters from session files into genbank or fasta files"""


import logging
import time
import requests

from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# make sure that pre and post 1.78 biopython are valid
try:
    from Bio.Alphabet import generic_dna
except ImportError:
    generic_dna = None
from Bio.SeqFeature import SeqFeature, FeatureLocation


from cblaster.classes import Session
from cblaster.extract import organism_matches, parse_organisms, parse_scaffolds
from cblaster.helpers import efetch_sequences
from cblaster.database import query_sequences, query_nucleotides


LOG = logging.getLogger(__name__)

# from https://www.ncbi.nlm.nih.gov/books/NBK25497/
MAX_REQUEST_SIZE = 500
MIN_TIME_BETWEEN_REQUEST = 0.34  # seconds


def parse_numbers(cluster_numbers):
    """Parses cluster numbers from user input.

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
            LOG.warning(f"Cannot extract cluster '{number}': number is not a valid integer")
    return chosen_cluster_numbers


def get_sorted_cluster_hierarchies(
    session,
    cluster_numbers=None,
    score_threshold=None,
    organisms=None,
    scaffolds=None,
    max_clusters=50,
):
    """Filter out selected clusters with their associated scaffold and organism

    Args:
        session (Session): A session object
        cluster_numbers (list): Numbers of clusters or a range of numbers eg 1-5
        score_threshold (float): Minimum score a cluster needs to have in order to be included
        organisms (list): Regex patterns for organisms of which all clusters need to be extracted
        scaffolds (list): Names of scaffolds of which all clusters need to be extracted
        max_clusters (int, None): the maximum amount of clusters extracted regardless of filters.
            You can set the value to None to extract all clusters
    Returns:
        List of tuples of in the form (cblaster.Cluster object, cblaster.Scaffold object of cluster,
         organism_name of cluster) sorted on cluster score
    """
    # No filter options return all clusters
    selected_clusters = set()

    # Prepare the filters defined by the user
    if cluster_numbers:
        cluster_numbers = parse_numbers(cluster_numbers)
    if organisms:
        organisms = parse_organisms(organisms)
    if scaffolds:
        scaffolds = parse_scaffolds(scaffolds)

    # Actually filter out the clusters
    for organism in session.organisms:
        if organisms and not organism_matches(organism.name, organisms):
            continue
        for scaffold in organism.scaffolds.values():
            if scaffolds and scaffold.accession not in scaffolds:
                continue
            for cluster in scaffold.clusters:
                if cluster_numbers and cluster.number not in cluster_numbers:
                    continue
                if score_threshold and cluster.score < score_threshold:
                    continue
                if scaffolds and not cluster_in_range(
                    scaffolds[scaffold.accession]["start"],
                    scaffolds[scaffold.accession]["end"],
                    cluster,
                ):
                    continue
                selected_clusters.add((cluster, scaffold, organism.full_name))

    # Make sure the sort is consistent and that same, scores, locations
    # are always sorted in the same way.
    return sorted(
        selected_clusters,
        key=lambda x: (x[0].score, -x[0].start, -x[0].end, x[1].accession),
        reverse=True,
    )[:max_clusters]


def cluster_in_range(start, end, cluster):
    """Checks if a cluster is within a given range.

    Args:
        start (int): start of the range
        end (int): end of the range
        cluster (Cluster): cblaster cluster object
    Returns:
        boolean given if the cluster is inbetween the start and end value
    """
    if start is end is None:
        return True
    return start <= cluster.start and end >= cluster.end


def create_genbanks_from_clusters(
    session, cluster_hierarchy, output_dir, prefix, format_
):
    """Create genbank files for each selected cluster

    Args:
        session (Session): a cblaster session object
        cluster_hierarchy (List): a sorted list of clusters with scaffold and organism
        output_dir (string): path to a directory for writing the output files
        prefix (string): string to start the file name of each cluster with
        format_ (str): the format that the extracted cluster should have
    """
    if session.params["mode"] == "remote":
        proteins = efetch_protein_sequences(cluster_hierarchy)
        nucleotides = efetch_nucleotide_sequence(cluster_hierarchy)
        name_attr = "name"
    elif session.params["mode"] == "local":
        sqlite_db = session.params["sqlite_db"]
        proteins = local_fetch_sequences(sqlite_db, cluster_hierarchy)
        nucleotides = local_fetch_nucleotide(sqlite_db, cluster_hierarchy)
        name_attr = "id"
    else:
        raise NotImplementedError(f"No protocol for mode {session.params['mode']}")

    output_dir = Path(output_dir)

    # Make the directory if it doesn't already exist
    if not output_dir.is_dir():
        output_dir.mkdir()

    # Generate genbank files for all the required clusters
    for cluster, scaffold, organism in cluster_hierarchy:
        cluster_proteins = {
            getattr(subject, name_attr): proteins[getattr(subject, name_attr)]
            for subject in [*cluster.subjects, *cluster.intermediate_genes]
        }
        output_file = output_dir / f"{prefix}cluster{cluster.number}.gbk"
        with output_file.open("w") as fp:
            record = cluster_to_record(
                cluster,
                cluster_proteins,
                nucleotides.get(cluster.number),
                organism,
                scaffold.accession,
                format_,
                session.params["require"],
            )
            SeqIO.write(record, fp, "genbank")

        LOG.debug(f"Created {output_file.name} file for cluster {cluster.number}")


def local_fetch_nucleotide(sqlite_db, cluster_hierarchy):
    """Fetches nucleotide sequence for clusters from SQLite3 database."""
    sequences = {}

    for cluster, scaffold, organism in cluster_hierarchy:
        sequence, *_ = query_nucleotides(
            scaffold.accession,
            organism,
            cluster.intermediate_start,
            cluster.intermediate_end,
            sqlite_db,
        )
        sequences[cluster.number] = sequence

    return sequences


def local_fetch_sequences(sqlite_db, cluster_hierarchy):
    """Fetches sequences from an offline SQLite3 database.

    Args:
        sqlite_db (str): path to a json_db file
        cluster_hierarchy (Set): set of tuples in the form (cblaster.Cluster object,
         scaffold_accession of cluster, organism_name of cluster)
    Returns:
        tuple of a dict of protein sequences keyed on protein_names and a dict of
         nucleotide sequences keyed on scaffold_accessions
    """
    subject_ids = [
        subject.id
        for cluster, _, _ in cluster_hierarchy
        for subject in [*cluster.subjects, *cluster.intermediate_genes]
    ]
    return {id: translation for id, translation in query_sequences(subject_ids, sqlite_db)}


def efetch_protein_sequences(cluster_hierarchy):
    """EFetches all protein sequences from all clusters in the cluster_hierarchy

    Args:
        cluster_hierarchy (Set): set of tuples in the form (cblaster.Cluster object,
         scaffold_accession of cluster, organism_name of cluster)
    Returns:
        a dictionary with protein sequences keyed on protein names
    """
    subject_names = [
        subject.name
        for cluster, _, _ in cluster_hierarchy
        for subject in [*cluster.subjects, *cluster.intermediate_genes]
    ]
    return efetch_sequences(subject_names)


def efetch_nucleotide_sequence(cluster_hierarchy):
    """Efetches all nucleotide sequences for all clusters in the cluster_hierarchy

    Args:
        cluster_hierarchy (Set):
            set of tuples in the form (cblaster.Cluster object,
            scaffold_accession of cluster, organism_name of cluster)
    Returns:
        a dictionary with nucleotide sequences keyed on cluster number
    """
    sequences = dict()
    passed_time = 0

    for cluster, scaffold, org_name in cluster_hierarchy:
        if passed_time < MIN_TIME_BETWEEN_REQUEST:
            time.sleep(MIN_TIME_BETWEEN_REQUEST - passed_time)
        start_time = time.time()

        LOG.info(
            f"Querying NCBI for sequence of {scaffold.accession}"
            f" from {cluster.intermediate_start}"
            f" to {cluster.intermediate_end}"
        )
        response = requests.post(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
            params={
                "db": "sequences",
                "rettype": "fasta",
                "seq_start": str(cluster.intermediate_start),
                "seq_stop": str(cluster.intermediate_end),
                "strand": "1",
            },
            files={"id": scaffold.accession},
        )
        LOG.debug(f"Efetch URL: {response.url}")

        if response.status_code != 200:
            raise requests.HTTPError(
                f"Error fetching sequences from NCBI [code {response.status_code}]."
                " Incorrect scaffold accession?"
            )

        # Only save the sequence not the fasta header
        sequences[cluster.number] = (
            response.text
            .strip()
            .split("\n", 1)[1]
            .replace("\n", "")
        )
        passed_time = time.time() - start_time

    return sequences


def cluster_to_record(
    cluster,
    cluster_prot_sequences,
    cluster_nuc_sequence,
    organism_name,
    scaffold_accession,
    format_,
    required_genes,
    mode="remote",
):
    """Convert a cblaster.Cluster object into a Bio.Seqrecord object

    Args:
        cluster (cblaster.Cluster): cblasdter Cluster object
        cluster_prot_sequences (dict): dictionary of protein sequences keyed on protein ids
        cluster_nuc_sequence (str): cluster nucleotide sequence
        organism_name (str): name of the organism the cluster is originated from
        scaffold_accession (str): accession of the scaffold the cluster belongs to
        format_ (str): the format that the extracted cluster should have
        required_genes (list): list of genes that the user has marked as required to be present in the cluster
    Returns:
        a Bio.Seqrecord object
    """
    if generic_dna:
        nuc_seq_obj = Seq(cluster_nuc_sequence, generic_dna)
    else:
        # Newer Biopython refuses second argument
        nuc_seq_obj = Seq(cluster_nuc_sequence)

    record = SeqRecord(
        nuc_seq_obj,
        id=scaffold_accession,
        name=scaffold_accession,
        annotations={"molecule_type": "DNA"},
        description=f"Genes for cluster {cluster.number} on scaffold {scaffold_accession} of species {organism_name}",
    )
    source_feature = SeqFeature(
        FeatureLocation(
            start=cluster.start - cluster.intermediate_start,
            end=cluster.end - cluster.intermediate_start),
        type="SOURCE",
        qualifiers={"organism": organism_name},
    )
    record.features.append(source_feature)

    if format_ == "bigscape":
        region_feature = SeqFeature(
            FeatureLocation(start=cluster.start - cluster.intermediate_start,
                            end=cluster.end - cluster.intermediate_start),
            type="region",
            qualifiers={"product": "other"},
        )
        record.features.append(region_feature)

    # Sequences in cluster_prot_sequences will be keyed on row ID if
    # they are retrieved from a local database, so use subject IDs if not null
    subjects = {
        subject.id if subject.id else subject.name: subject
        for subject in [*cluster.subjects, *cluster.intermediate_genes]
    }
    for key, sequence in cluster_prot_sequences.items():
        subject = subjects[key]
        qualifiers = {"protein_id": subject.name, "translation": sequence}

        top_hit = None
        if len(subject.hits) > 0:
            top_hit = max(subject.hits, key=lambda x: x.bitscore)

        # Indicate the role of a gene in a cluster:
        #   hit_required = hit with a required gene
        #   hit = a hit with a non required gene
        #   intermediate = an intermediate gene
        if top_hit is not None:
            if required_genes is None or top_hit.query in required_genes:
                qualifiers["cluster_role"] = "hit_required"
                if format_ == "bigscape":
                    qualifiers["gene_kind"] = "biosynthetic"
            else:
                qualifiers["cluster_role"] = "hit"
        else:
            qualifiers["cluster_role"] = "intermediate"

        # Build the SeqFeature object corresponding to the CDS
        location = FeatureLocation(
            start=subject.start - cluster.intermediate_start,
            end=subject.end - cluster.intermediate_start,
            strand=1 if subject.strand == "+" else -1,
        )
        cds_feature = SeqFeature(location, type="CDS", qualifiers=qualifiers)
        record.features.append(cds_feature)

    record.features.sort(key=lambda x: x.location.start)
    return record


def extract_clusters(
    session,
    output_dir,
    prefix="",
    cluster_numbers=None,
    score_threshold=None,
    organisms=None,
    scaffolds=None,
    format_="genbank",
    max_clusters=50,
):
    """Extracts Cluster objects from a Session file and writes them to a file.

    If BiG-SCAPE format is chosen,  a 'gene_kind' qualifier is added to each CDS
    feature to indicate what genes are part of the core of the cluster.

    Genes that are flagged as required are considered core. If no genes are flagged as
    required, all genes are considered to be core genes.

    An additional qualifier is provided called 'cluster_role'. This qualifier allows
    the identification between hits found against required genes of the query, hits
    found agains any gene of the query and intermediate genes.

    Args:
        session (string): path to a session.json file
        output_dir (string): path to a directory for writing the output files
        prefix (string): string to start the file name of each cluster with
        cluster_numbers (list): cluster numbers to include
        score_threshold (float): minum score in order for a cluster to be included
        organisms (list): Organism filtering regular expressions, clusters for
         these organisms are included
        scaffolds(list): clusters on these scaffolds are included
        format_ (str): the format that the extracted cluster should have
        max_clusters (int): the maximum amount of clusters extracted regardless of filters

    """
    LOG.info("Starting cblaster plotting of clusters using clinker")
    LOG.info("Loading session from: %s", session)
    with open(session) as fp:
        session = Session.from_json(fp)

    LOG.info("Extracting clusters that match the filters")
    cluster_hierarchy = get_sorted_cluster_hierarchies(
        session,
        cluster_numbers,
        score_threshold,
        organisms,
        scaffolds,
        max_clusters,
    )

    LOG.info(f"Extracted {len(cluster_hierarchy)} clusters.")
    if len(cluster_hierarchy) == 0:
        LOG.info("There are no clusters that meet the filtering criteria. Exiting...")
        raise SystemExit(0)

    LOG.info("Writing genbank files")
    create_genbanks_from_clusters(
        session,
        cluster_hierarchy,
        output_dir,
        prefix,
        format_,
    )

    LOG.info(f"Clusters have been written to {output_dir}")
    LOG.info("Done!")
