#!/usr/bin/env python3

"""Add intermediate genes to the clusters of a session"""


import logging
import time
import requests
import re
import io

from Bio import SeqIO

from cblaster.extract_clusters import get_sorted_cluster_hierarchies
from cblaster.database import query_intermediate_genes
from cblaster.classes import Subject
from cblaster.genome_parsers import seqrecord_to_tuples


LOG = logging.getLogger(__name__)

# from https://www.ncbi.nlm.nih.gov/books/NBK25497/
MIN_TIME_BETWEEN_REQUEST = 0.34  # seconds
PROTEIN_NAME_IDENTIFIERS = ("protein_id", "locus_tag", "gene", "ID", "Name", "label")


def set_local_intermediate_genes(sqlite_db, cluster_hierarchy, gene_distance):
    """Adds intermediate genes to clusters in the cluster_hierarchy using a SQLite database

    Args:
        sqlite_db (str): path to the sqlite database
        cluster_hierarchy (List): Tuples with cblaster cluster scaffold accession and organism name
        gene_distance (int): the extra distance around a cluster to collect genes from
    """
    for cluster, scaffold, organism in cluster_hierarchy:
        scaffold_accession = scaffold.accession
        search_start = cluster.start - gene_distance
        search_stop = cluster.end + gene_distance
        cluster_ids = [subject.id for subject in cluster.subjects]
        cluster.intermediate_genes = [
            Subject(id=id, name=name, start=start, end=end, strand=strand)
            for start, end, id, name, strand in query_intermediate_genes(
                cluster_ids,
                search_start,
                search_stop,
                scaffold_accession,
                organism,
                sqlite_db,
            )
        ]


def intermediate_genes_request(accession, start, stop):
    response = requests.post(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
        params={
            "db": "nuccore",
            "rettype": "gb",
            "from": str(start),
            "to": str(stop),
        },
        files={"id": accession},
    )
    LOG.info(f"Fetching intermediate genes from NCBI from {accession}")
    LOG.debug(f"Efetch URL: {response.url}")
    if response.status_code != 200:
        raise requests.HTTPError(
            f"Error fetching intermediate genes for NCBI [code {response.status_code}]."
        )
    return response


def intermediate_genes_genbank(response):
    features = []
    for record in SeqIO.parse(io.StringIO(response.text), "genbank"):
        subjects = []
        for feature in seqrecord_to_tuples(record, None):
            try:
                subject = tuple_to_subject(feature)
            except TypeError:
                continue
            subjects.append(subject)
        features.extend(subjects)
    return features


def tuple_to_subject(feature):
    assert len(feature) == 8, "Tuple should be of length 8"
    ftype, name, start, end, strand, sequence, record_id, source_id = feature
    if ftype != "gene":
        raise TypeError(f"Feature type {ftype}, expecting 'gene'")
    return Subject(
        name=name,
        start=start,
        end=end,
        strand=strand,
        sequence=sequence,
    )


def set_remote_intermediate_genes(cluster_hierarchy, gene_distance):
    """Adds intermediate genes to clusters in the cluster_hierarchy from NCBI feature tables.

    Args:
        cluster_hierarchy (List): Tuples with Cblaster cluster scaffold accession and organism name
        gene_distance (int): the extra distance around a cluster to collect genes from
    """
    passed_time = 0
    for cluster, scaffold, _ in cluster_hierarchy:
        scaffold_accession = scaffold.accession
        if passed_time < MIN_TIME_BETWEEN_REQUEST:
            time.sleep(MIN_TIME_BETWEEN_REQUEST - passed_time)
        search_start = max(0, cluster.start - gene_distance)
        search_stop = cluster.end + gene_distance
        start_time = time.time()
        response = intermediate_genes_request(scaffold_accession, search_start, search_stop)
        subjects = intermediate_genes_genbank(response)
        cluster.intermediate_genes = get_remote_intermediate_genes(subjects, cluster)
        passed_time = time.time() - start_time


def get_remote_intermediate_genes(subjects, cluster):
    """
    Get all genes that are not part of the cluster from the list of subjects

    Args:
        subjects (List): list of cblaster Subject objects
        cluster (Cluster): cblaster Cluster object

    Returns:
        List of all the genes that are in subjects and not in the cluster
    """
    cluster_genes = set([s.name for s in cluster.subjects])
    intermediate_genes = []
    for subject in subjects:
        if subject.name not in cluster_genes:
            intermediate_genes.append(subject)
    return intermediate_genes


def find_intermediate_genes(session, gene_distance=5000, max_clusters=100):
    """
    Main function called for finding intermediate genes.

    Args:
        session (Session): cblaster Session object
        gene_distance (int): the extra distance around a cluster to collect genes from
        max_clusters (int): maximum amount of clusters intermediate genes are added to
        considering that retrieving intermediate genes for remote sessions can become
        expensive.
    """
    LOG.info("Searching for intermediate genes")
    cluster_hierarchy = get_sorted_cluster_hierarchies(session, max_clusters=max_clusters)

    if session.params["mode"] == "local":
        set_local_intermediate_genes(
            session.params["sqlite_db"], cluster_hierarchy, gene_distance
        )
    elif session.params["mode"] == "remote":
        set_remote_intermediate_genes(cluster_hierarchy, gene_distance)
    else:
        LOG.warning(
            f"{session.params['mode']} is not supported for intermediated genes."
            f" Skipping intermediate genes addition"
        )
