
import logging

from cblaster.extract_clusters import extract_cluster_hierarchies
from cblaster.database import query_database_for_intermediate_genes


LOG = logging.getLogger(__name__)
MAX_CLUSTER_DISTANCE = 5000
PROTEIN_NAME_IDENTIFIERS = ("protein_id", "locus_tag", "gene", "ID", "Name", "label")


def get_local_intermediate_hits(sqlite_db, cluster_hierarchy):
    for cluster, _, _ in cluster_hierarchy:
        search_start, search_stop = cluster.start - MAX_CLUSTER_DISTANCE, cluster.end + MAX_CLUSTER_DISTANCE
        cluster_ids = [subject.name for subject in cluster.subjects]

        intermediate_genes = []
        for start, end, name, strand in \
                query_database_for_intermediate_genes(cluster_ids, search_start, search_stop, sqlite_db):
            intermediate_genes.append({
                "start": start,
                "end": end,
                "name": name,
                "strand": strand
            })
        cluster.intermediate_genes = intermediate_genes


def find_intermediate_hits(session):
    LOG.info("Searching for intermediate genes")
    cluster_hierarchy = extract_cluster_hierarchies(session)

    if session.params["mode"] == "local":
        sqlite_db = session.params["sqlite_db"]
        get_local_intermediate_hits(sqlite_db, cluster_hierarchy)
    elif session.params["mode"] == "remote":
        pass



