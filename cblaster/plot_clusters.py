
from pathlib import Path
import logging
from g2j import genbank
from clinker.classes import (
    Cluster as ClinkerCluster,
    Locus as ClinkerLocus,
    Gene as ClinkerGene
)
from clinker.align import (
    Alignment as ClinkerAlignment,
    Globaligner as ClinkerGlobalaligner
)
from clinker.plot import plot_clusters as clinker_plot_clusters


from cblaster.extract_clusters import extract_cluster_hierarchies
from cblaster.classes import Session
from cblaster import embl


LOG = logging.getLogger(__name__)


def query_to_clinker_cluster(query_file):
    with open(query_file) as query:
        if any(query_file.endswith(ext) for ext in (".gbk", ".gb", ".genbank", ".gbff")):
            organism = genbank.parse(query, feature_types=["CDS"])
        elif any(query_file.endswith(ext) for ext in (".embl", ".emb")):
            organism = embl.parse(query_file, feature_types=["CDS"])
        # TODO add the fasta case

    identifiers = ("protein_id", "locus_tag", "gene", "ID", "Name", "label")

    loci = []
    count = 1
    for locus_nr, scaffold in enumerate(organism.scaffolds):
        locus_genes = []
        sorted_cds_features = sorted(scaffold.features, key=lambda f: f.location.min())
        for feature in sorted_cds_features:
            name = None
            for identifier in identifiers:
                if identifier in feature.qualifiers:
                    name = feature.qualifiers[identifier].split(" ")[0]
                    break
            if not name:
                name = f"protein_{count}"
                count += 1

            locus_genes.append(ClinkerGene(label=name, start=feature.location.min(), end=feature.location.max(),
                                           strand=1 if feature.location.strand == '+' else -1))
        loci.append(ClinkerLocus(f"Locus{locus_nr}", locus_genes, start=sorted_cds_features[0].location.min(),
                                 end=sorted_cds_features[-1].location.max()))
    return ClinkerCluster("Query_cluster", loci)


def clusters_to_clinker_alignments(clinker_query_cluster, clusters):
    """Create clinker.Alignments classes between the query cluster and all other clusters

    Make clinker.Link objects between all genes of the query and the genes in the clusters that
    where matched during blasting.

    Args:
        clinker_query_cluster(clinker.Cluster): clinker.Cluster object of the query used for the session
        clusters(List): a list of cblaster.Cluster objects
    Returns:
        a list of clinker.Alignment objects
    """
    allignments = []
    for cblaster_cluster in clusters:
        clinker_cluster = cblaster_cluster.to_clinker_cluster()
        allignment = ClinkerAlignment(query=clinker_query_cluster, target=clinker_cluster)
        for subject in cblaster_cluster.subjects:
            best_hit = max(subject.hits, key=lambda x: x.bitscore)
            query_gene = _gene_from_clinker_cluster(clinker_query_cluster, best_hit.query)
            subject_gene = _gene_from_clinker_cluster(clinker_cluster, best_hit.subject)
            allignment.add_link(query_gene, subject_gene, best_hit.identity, 0)
        allignments.append(allignment)
    return allignments


def _gene_from_clinker_cluster(cluster, gene_label):
    """Get a gene from a clinker.Cluster object

    works the same as clinker.Cluster.get_gene, except that this function is
    broken at the moment since clinker.Gene objects have no name attribute
    anymore but label objects instead.

    Args:
        cluster(clinker.Cluster): a clinker.Cluster object
        gene_label(str): the label of the clinker.Gene object that is requested
    Returns:
        a clinker.Gene object or None of no sutch gene exists in the cluster.
    """
    for locus in cluster.loci:
        for gene in locus.genes:
            if gene.label == gene_label:
                return gene


def allignments_to_clinker_global_alligner(allignments):
    """Create a clinker.Globalaligner object from alignments

    Args:
        allignments (List): a list of clinker.Aligner objects
    Returns:
        a clinker.Globalaligner object
    """
    global_aligner = ClinkerGlobalaligner()
    for allignment in allignments:
        global_aligner.add_alignment(allignment)
    return global_aligner


def plot_clusters(
    session,
    cluster_numbers=None,
    score_threshold=None,
    organisms=None,
    scaffolds=None,
    plot_outfile=None,
):
    """Plot Cluster objects from a Session file

    Args:
        session (string): path to a session.json file
        cluster_numbers (list): cluster numbers to include
        score_threshold (float): minum score in order for a cluster to be included
        organisms (list): Organism filtering regular expressions, clusters for
         these organisms are included
        scaffolds(list): clusters on these scaffolds are included
        plot_outfile (str): path to a file for the final plot
    """
    logging.info("Starting generation of cluster plot with clinker.")
    with open(session, "r") as f:
        session = Session.from_json(f.read())

    # filter the cluster using the filter functions from extract_clusters modue
    cluster_hierarchies = extract_cluster_hierarchies(session, cluster_numbers, score_threshold, organisms, scaffolds)

    clinker_query_cluster = query_to_clinker_cluster(session.params["query_file"])

    allignments = clusters_to_clinker_alignments(clinker_query_cluster, [h[0] for h in cluster_hierarchies])
    global_aligner = allignments_to_clinker_global_alligner(allignments)

    clinker_plot_clusters(global_aligner, plot_outfile, use_file_order=True)
    if plot_outfile:
        LOG.info(f"Plot file can be found at {plot_outfile}")
    LOG.info("Done!")
