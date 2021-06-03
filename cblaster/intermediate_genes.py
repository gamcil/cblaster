#!/usr/bin/env python3

"""Add intermediate genes to the clusters of a session"""


import logging
import time
import requests
import re

from cblaster.extract_clusters import get_sorted_cluster_hierarchies
from cblaster.database import query_intermediate_genes
from cblaster.classes import Subject


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
            "rettype": "ft",
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
        response = intermediate_genes_request(
            scaffold_accession, search_start, search_stop
        )
        subjects = genes_from_feature_table(response.text, search_start)
        cluster.intermediate_genes = get_remote_intermediate_genes(subjects, cluster)
        passed_time = time.time() - start_time


def transfer_protein_ids(genes, cdss):
    """Transfers protein_id qualifier from CDS to matching gene features."""
    for gene in genes:
        for cds in cdss:
            if gene["start"] <= cds["start"] and cds["end"] <= gene["end"]:
                gene["protein_id"] = cds["protein_id"]
                break


def parse_table_features(table_text):
    """Parse features into dictionaries containing start, end and qualifiers """
    features = []
    feature = {}
    for line in table_text.split("\n"):
        tabs = line.split("\t")
        if len(tabs) == 3:
            if feature:
                features.append(feature)
            start, end, ftype = tabs
            start, end, strand = get_start_end_strand(start, end)
            feature = dict(
                type=ftype,
                start=start,
                end=end,
                strand=strand,
                locus_tag=None,
                protein_id=None,
            )
        elif len(tabs) == 2 and feature:
            start, end = tabs
            start, end, _ = get_start_end_strand(start, end) 
            feature["start"] = min(feature["start"], start)
            feature["end"] = max(feature["end"], end)
        elif len(tabs) == 5:
            *_, key, value = tabs
            if key in ("locus_tag", "protein_id"):
                feature[key] = value

    # Make sure we catch the last feature
    if feature:
        features.append(feature)

    # If there are gene features, prefer using their coordinates instead of CDS.
    # Just transfer protein_ids from matching CDS features and return.
    # Otherwise, just return CDS features.
    genes = [f for f in features if f["type"]  == "gene"]
    cdss = [f for f in features if f["type"]  == "CDS"]
    if genes:
        transfer_protein_ids(genes, cdss)
        return genes
    return cdss


def set_gene_name(feature):
    """Sets the name of a gene feature.

    Tries to extract protein_id from feature table protein_id qualifier,
    which takes the form:
        source database|protein_id|additional qualifiers

    For example:
        ref|WP_063781455.1|
        gb|ACG70812.1|
        emb|SPB51992.1||gnl|WGS:OGUI|SPB51992
        dbj|GAQ43934.1|

    If this fails, falls back to locus_tag.
    """
    pattern = re.compile(r"^\w+\|(?P<protein_id>[A-Za-z0-9\._-]+)(?:\||$)")
    protein_id = feature.pop("protein_id", None)
    locus_tag = feature.pop("locus_tag", None)
    try:
        feature["name"] = pattern.search(protein_id).group("protein_id")
    except:
        feature["name"] = locus_tag


def genes_from_feature_table(table_text, search_start):
    """Extracts all CDS regions from a feature table returned by NCBI.

    Additional information for the feature table can be found here:
    https://www.ncbi.nlm.nih.gov/WebSub/html/help/feature-table.html

    Args:
        table_text (str): The feature table in text format given by NCBI
        search_start (int): The base pair start of the region where genes are returned from
            since the location of genes provided is relative to the requested region.
    Returns:
        List of cblaster Subject objects containing the intermediate genes.
    """

    # Parse features into dictionaries containing start, end and qualifiers
    features = parse_table_features(table_text)

    # Instantiate Subject objects
    # - Add scaffold search start to gene start/end
    # - Extract protein IDs from protein_id qualifiers
    subjects = []
    for feature in features:
        feature.pop("type")
        feature["start"] += search_start
        feature["end"] += search_start
        set_gene_name(feature)
        subject = Subject(**feature)
        subjects.append(subject)

    return subjects


def get_start_end_strand(start, end):
    """Extracts the start and end locations from the start and end string.

    The location can contain < and > which should be removed. Additionally
    if the end is smaller then the start the location is on the negative strand

    Args:
        start(str): string representing the start location of a gene
        end(str): string representing the end location of a gene

    Returns:
        a tuple with start as integer, end as integer and strand as '+' or '-'
    """
    start = int(start.replace("<", "").replace(">", ""))
    end = int(end.replace("<", "").replace(">", ""))
    strand = 1
    if start > end:
        strand = -1
        new_start = end
        end = start
        start = new_start
    return start, end, strand


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
