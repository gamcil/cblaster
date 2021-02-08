
import logging
import time
import requests
import re

from Bio import SeqFeature

from cblaster.extract_clusters import extract_cluster_hierarchies
from cblaster.database import query_database_for_intermediate_genes
from cblaster.classes import Subject


LOG = logging.getLogger(__name__)
PROTEIN_NAME_IDENTIFIERS = ("protein_id", "locus_tag", "gene", "ID", "Name", "label")

# from https://www.ncbi.nlm.nih.gov/books/NBK25497/
MIN_TIME_BETWEEN_REQUEST = 0.34  # seconds


def set_local_intermediate_genes(sqlite_db, cluster_hierarchy, gene_distance):
    for cluster, _, _ in cluster_hierarchy:
        search_start, search_stop = cluster.start - gene_distance, cluster.end + gene_distance
        cluster_ids = [subject.name for subject in cluster.subjects]

        intermediate_genes = []
        for start, end, name, strand in \
                query_database_for_intermediate_genes(cluster_ids, search_start, search_stop, sqlite_db):
            # generate an empty subject
            intermediate_genes.append(Subject(name=name, start=start, end=end, strand="+" if strand == 1 else "-"))
        cluster.intermediate_genes = intermediate_genes


def set_remote_intermediate_genes(cluster_hierarchy, gene_distance):
    sequences = dict()
    passed_time = 0
    for cluster, scaffold_accession, _ in cluster_hierarchy:
        if passed_time < MIN_TIME_BETWEEN_REQUEST:
            time.sleep(MIN_TIME_BETWEEN_REQUEST - passed_time)
        search_start, search_stop = cluster.start - gene_distance, cluster.end + gene_distance
        start_time = time.time()
        response = requests.post(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
            params={"db": "nuccore", "rettype": "ft", "from": str(search_start), "to": str(search_stop)},
            files={"id": scaffold_accession}
        )
        LOG.info(f"Fetching intermediate genes from NCBI from {scaffold_accession}")
        LOG.debug(f"Efetch URL: {response.url}")
        if response.status_code != 200:
            raise requests.HTTPError(
                f"Error fetching intermediate genes from NCBI [code {response.status_code}]."
            )

        subjects = genes_from_feature_table(response.text)
        intermediate_genes = get_intermediate_genes(subjects, cluster)
        cluster.intermediate_genes = intermediate_genes
        passed_time = time.time() - start_time
    return sequences


def genes_from_feature_table(table_text):
    start = end = name = strand = None
    intermediate_genes = []
    for line in table_text.split("\n"):
        tabs = line.split("\t")
        if len(tabs) < 3:
            continue
        elif tabs[2] == "CDS":
            if name is not None:
                intermediate_genes.append(Subject(name=name, start=start, end=end, strand=strand))
                name = None
            start, end, strand = get_start_end_strand(tabs[0], tabs[1])
        elif len(tabs) == 5 and tabs[3] == "protein_id":
            name = re.search(r"\|([A-Za-z0-9\._]+)\|", tabs[4]).group(1)

    return intermediate_genes

#https://www.ncbi.nlm.nih.gov/WebSub/html/help/feature-table.html
def get_start_end_strand(start, end):
    start = int(start.replace("<", "").replace(">", ""))
    end = int(end.replace("<", "").replace(">", ""))
    strand = "+"
    if start > end:
        strand = "-"
        new_start = end
        end = start
        start = new_start
    return start, end, strand


def get_intermediate_genes(subjects, cluster):
    cluster_genes = set([s.name for s in cluster.subjects])
    intermediate_genes = []
    for subject in subjects:
        if subject.name not in cluster_genes:
            intermediate_genes.append(subject)
    return intermediate_genes


def find_intermediate_genes(session, gene_distance=5000, max_clusters=100):
    LOG.info("Searching for intermediate genes")

    cluster_hierarchy = extract_cluster_hierarchies(session, max_clusters=max_clusters)

    if session.params["mode"] == "local":
        sqlite_db = session.params["sqlite_db"]
        set_local_intermediate_genes(sqlite_db, cluster_hierarchy, gene_distance)
    elif session.params["mode"] == "remote":
        set_remote_intermediate_genes(cluster_hierarchy, gene_distance)
    else:
        LOG.warning(f"{session.params['mode']} is not supported for intermediated genes.")
