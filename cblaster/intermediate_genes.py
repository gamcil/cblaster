
import logging
import json
from g2j.classes import Organism

from cblaster.extract_clusters import extract_cluster_hierarchies


LOG = logging.getLogger(__name__)
MAX_CLUSTER_DISTANCE = 5000
PROTEIN_NAME_IDENTIFIERS = ("protein_id", "locus_tag", "gene", "ID", "Name", "label")


def get_intermediate_gene(feature, cluster_subject_names):
    name = None
    for identifier in PROTEIN_NAME_IDENTIFIERS:
        if identifier not in feature.qualifiers:
            continue
        elif feature.qualifiers[identifier] in cluster_subject_names:
            return None
        elif name is None:
            name = feature.qualifiers[identifier]
    if name is not None:
        return {
            "start": feature.location.min(),
            "end": feature.location.max(),
            "name": name,
            "strand": feature.location.strand
        }
    return None


def get_local_intermediate_hits(g2j_organisms, cluster_hierarchy):
    organisms_dict = {organism.name: organism for organism in g2j_organisms}

    for cluster, scaffold_accession, organism_name in cluster_hierarchy:
        cluster_boundaries = cluster.start - MAX_CLUSTER_DISTANCE, cluster.end + MAX_CLUSTER_DISTANCE
        cluster_subject_names = set(subject.name for subject in cluster.subjects)
        organism = organisms_dict[organism_name]
        intermediate_genes = []
        for scaffold in organism.scaffolds:
            if scaffold.accession != scaffold_accession:
                continue
            for feature in scaffold.features:
                if feature.type != "CDS" or \
                        feature.location.min() < cluster_boundaries[0] or \
                        feature.location.max() > cluster_boundaries[1]:
                    continue
                intermediate_gene = get_intermediate_gene(feature, cluster_subject_names)
                if intermediate_gene is not None:
                    intermediate_genes.append(intermediate_gene)
        cluster.intermediate_genes = intermediate_genes


def find_intermediate_hits(session):
    json_db = session.params["json_db"]
    # read the json database
    with open(json_db, "r") as f:
        # g2j Organism objects not cblaster Organism objects
        organism_objs = [Organism.from_dict(dct) for dct in json.load(f)]
    cluster_hierarchy = extract_cluster_hierarchies(session)

    if session.params["mode"] == "local":
        get_local_intermediate_hits(organism_objs, cluster_hierarchy)
    elif session.params["mode"] == "remote":
        pass



