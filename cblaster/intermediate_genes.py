
import logging

from cblaster.extract_clusters import extract_cluster_hierarchies
from cblaster.database import query_database_for_intermediate_genes
from cblaster.classes import Subject


LOG = logging.getLogger(__name__)
PROTEIN_NAME_IDENTIFIERS = ("protein_id", "locus_tag", "gene", "ID", "Name", "label")


def get_local_intermediate_genes(sqlite_db, cluster_hierarchy, gene_distance):
    for cluster, _, _ in cluster_hierarchy:
        search_start, search_stop = cluster.start - gene_distance, cluster.end + gene_distance
        cluster_ids = [subject.name for subject in cluster.subjects]

        intermediate_genes = []
        for start, end, name, strand in \
                query_database_for_intermediate_genes(cluster_ids, search_start, search_stop, sqlite_db):
            # generate an empty subject
            intermediate_genes.append(Subject(name=name, start=start, end=end, strand="+" if strand == 1 else "-"))
        cluster.intermediate_genes = intermediate_genes


def find_intermediate_genes(session, gene_distance=5000, max_clusters=100):
    LOG.info("Searching for intermediate genes")
    cluster_hierarchy = extract_cluster_hierarchies(session, max_clusters=max_clusters)

    if session.params["mode"] == "local":
        sqlite_db = session.params["sqlite_db"]
        get_local_intermediate_genes(sqlite_db, cluster_hierarchy, gene_distance)
    elif session.params["mode"] == "remote":
        pass
