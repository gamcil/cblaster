#!/usr/bin/env python3

"""Plot clusters from session files"""


import logging
from Bio import SeqIO

from clinker.classes import (
    Cluster as ClinkerCluster,
    Locus as ClinkerLocus,
    Gene as ClinkerGene
)
from clinker.align import (
    Alignment as ClinkerAlignment,
    Globaligner as ClinkerGlobalaligner,
    Group as ClinkerGroup,
)
from clinker.plot import plot_clusters as clinker_plot_clusters
from cblaster.extract_clusters import get_sorted_cluster_hierarchies
from cblaster.classes import Session
from cblaster.genome_parsers import GBK_SUFFIXES, EMBL_SUFFIXES


LOG = logging.getLogger(__name__)

FASTA_SPACE = 500


def cblaster_to_clinker_cluster(
    cluster,
    cluster_label=None,
    scaffold_accession="",
    organism_name="",
):
    """Convert this cluster to a clinker format cluster

    Args:
        scaffold_accession (str): accession of the scaffold this cluster is located on
    Returns:
        A clinker.Cluster object
    """
    clinker_genes = []
    for subject in cluster.subjects:
        tooltip_dict = dict(accession=subject.name)
        if subject.hits:
            best_hit = max(subject.hits, key=lambda x: x.bitscore)
            for key in ["identity", "coverage", "evalue", "bitscore"]:
                value = getattr(best_hit, key)
                if value:
                    tooltip_dict[key] = value
        clinker_gene = ClinkerGene(
            label=subject.name,
            start=subject.start,
            end=subject.end,
            strand=subject.strand,
            names=tooltip_dict
        )
        clinker_genes.append(clinker_gene)

    for gene in cluster.intermediate_genes:
        tooltip_dict = {"accession": gene.name}
        clinker_gene = ClinkerGene(
            label=gene.name,
            start=gene.start,
            end=gene.end,
            strand=gene.strand,
            names=tooltip_dict
        )
        clinker_genes.append(clinker_gene)

    clinker_locus = ClinkerLocus(
        scaffold_accession,
        clinker_genes,
        start=cluster.intermediate_start,
        end=cluster.intermediate_end
    )

    if cluster_label:
        label = cluster_label
    else:
        inner = (
            f"Cluster {cluster.number}" + f", {cluster.score:.2f} score"
            if cluster.score
            else ""
        )
        label = f"{organism_name} ({inner})"

    return ClinkerCluster(label, [clinker_locus])


def clusters_to_clinker_globaligner(clinker_query_cluster, cluster_hierarchies):
    """Create clinker.Alignments classes between the query cluster and all other clusters

    Make clinker.Link objects between all genes of the query and the genes in the clusters that
    where matched during blasting.

    Args:
        clinker_query_cluster (clinker.Cluster):
            clinker.Cluster object of the query used for the session
        cluster_hierarchies (List):
            a list of tuples in the form (cblaster.Cluster object, scaffold_accession
            of cluster, organism_name of cluster)
    Returns:
        a list of clinker.Alignment objects
    """
    globaligner = ClinkerGlobalaligner()

    # Form gene groupings based on query genes
    groups = {
        gene.label: ClinkerGroup(label=gene.label, genes=[gene.uid])
        for locus in clinker_query_cluster.loci
        for gene in locus.genes
    }

    # Iterate all cblaster Clusters, convert to clinker Clusters and mock alignments
    for cblaster_cluster, scaffold, organism_name in cluster_hierarchies:
        clinker_cluster = cblaster_to_clinker_cluster(
            cblaster_cluster,
            scaffold_accession=scaffold.accession,
            organism_name=organism_name,
        )
        alignment = ClinkerAlignment(query=clinker_query_cluster, target=clinker_cluster)

        for subject in cblaster_cluster.subjects:
            # Find the best query hit
            best_hit = max(subject.hits, key=lambda x: x.bitscore)

            # Pull out the corresponding query/subject genes and create alignment
            query_gene = clinker_query_cluster.get_gene(best_hit.query)
            subject_gene = clinker_cluster.get_gene(best_hit.subject)
            alignment.add_link(query_gene, subject_gene, best_hit.identity / 100, 0)

            # Save the UID of the subject to the corresponding query group
            groups[query_gene.label].genes.append(subject_gene.uid)

        globaligner.add_alignment(alignment)
    globaligner.groups = [group for group in groups.values()]
    return globaligner


def plot_clusters(
    session,
    cluster_numbers=None,
    score_threshold=None,
    organisms=None,
    scaffolds=None,
    plot_outfile=None,
    max_clusters=50,
    testing=False,
):
    """Plot Cluster objects from a Session file
    Args:
        session (string): path to a session.json file
        cluster_numbers (list): cluster numbers to include
        score_threshold (float): minumum score in order for a cluster to be included
        organisms (list): Organism filtering regular expressions, clusters for
        these organisms are included
        scaffolds(list): clusters on these scaffolds are included
        plot_outfile (str): path to a file for the final plot
        max_clusters (int): the maximum amount of clusters plotted regardless of filters
        testing (bool): argument to switch of plotting when testing making sure that no dynamioc plot
        is served since this will crash the testing.
    """
    LOG.info("Starting generation of cluster plot with clinker.")
    session = Session.from_file(session)

    # Filter the cluster using filter functions from the extract_clusters module
    cluster_hierarchies = get_sorted_cluster_hierarchies(
        session,
        cluster_numbers,
        score_threshold,
        organisms,
        scaffolds,
        max_clusters,
    )

    # Form the query cluster from the session query file
    query_cluster = cblaster_to_clinker_cluster(
        session.query,
        cluster_label="Query Cluster",
        scaffold_accession=session.params.get("query_file", "N.A."),
    )

    # Create a Globaligner object containing mocked clusters/alignments/links
    globaligner = clusters_to_clinker_globaligner(query_cluster, cluster_hierarchies)

    if not testing:
        clinker_plot_clusters(globaligner, plot_outfile, use_file_order=True)

    if plot_outfile:
        LOG.info(f"Plot file can be found at {plot_outfile}")

    LOG.info("Done!")
