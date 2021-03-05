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


def cblaster_to_clinker_cluster(cluster, scaffold_accession="", organism_name=""):
    """Convert this cluster to a clinker format cluster

    Args:
        scaffold_accession (str): accession of the scaffold this cluster is located on
    Returns:
        A clinker.Cluster object
    """
    clinker_genes = []
    for subject in cluster.subjects:
        best_hit = max(subject.hits, key=lambda x: x.bitscore)
        tooltip_dict = {
            "accession": subject.name,
            "identity": f"{best_hit.identity:.2f}",
            "bitscore": best_hit.bitscore,
            "coverage": f"{best_hit.coverage:.2f}",
            "e-value": best_hit.evalue if best_hit.evalue != 0 else "0.0"
        }
        clinker_gene = ClinkerGene(
            label=subject.name,
            start=subject.start,
            end=subject.end,
            strand=1 if subject.strand == '+' else -1,
            names=tooltip_dict
        )
        clinker_genes.append(clinker_gene)

    for gene in cluster.intermediate_genes:
        tooltip_dict = {"accession": gene.name}
        clinker_gene = ClinkerGene(
            label=gene.name,
            start=gene.start,
            end=gene.end,
            strand=1 if gene.strand == '+' else -1,
            names=tooltip_dict
        )
        clinker_genes.append(clinker_gene)

    clinker_locus = ClinkerLocus(
        scaffold_accession,
        clinker_genes,
        start=cluster.intermediate_start,
        end=cluster.intermediate_end
    )

    return ClinkerCluster(
        f"{organism_name} Cluster {cluster.number} ({cluster.score:.2f} score)",
        [clinker_locus],
    )


def query_to_clinker_cluster(query_file):
    """Turn a query file of a cblaster Session object into a clinker.Cluster object
    Args:
        query_file (str): Path to the query file
    Returns:
        a clinker.Cluster object
    """
    with open(query_file) as query:
        if any(query_file.endswith(ext) for ext in GBK_SUFFIXES):
            seqrecord = SeqIO.parse(query, "genbank")
        elif any(query_file.endswith(ext) for ext in EMBL_SUFFIXES):
            seqrecord = SeqIO.parse(query, "embl")
        else:
            return fasta_to_cluster(query)
        return _seqrecord_to_clinker_cluster(seqrecord)


def _seqrecord_to_clinker_cluster(seqrecord):
    """Transform a Bio.Seqrecord into a clinker cluster so clinker can plot it easily
    Args:
        seqrecord (Bio.Seqrecord): seqrecord object
    Returns:
        A cinker cluster object
    """
    identifiers = ("protein_id", "locus_tag", "gene", "ID", "Name", "label")
    loci = []
    count = 1
    for locus_nr, record in enumerate(seqrecord):
        locus_genes = []
        locus_start = locus_end = None
        for feature in record.features:
            if feature.type != "CDS":
                continue
            name = None
            for identifier in identifiers:
                if identifier not in feature.qualifiers:
                    continue
                name = feature.qualifiers[identifier][0]
                break
            if name is None:
                name = f"protein_{count}"
                count += 1
            if locus_start is None or feature.location.start < locus_start:
                locus_start = feature.location.start
            if locus_end is None or feature.location.end > locus_end:
                locus_end = feature.location.end
            gene = ClinkerGene(
                label=name,
                start=feature.location.start,
                end=feature.location.end,
                strand=feature.location.strand,
                names={"accession": name}
            )
            locus_genes.append(gene)
        locus = ClinkerLocus(f"Locus{locus_nr + 1}", locus_genes, start=locus_start, end=locus_end)
        loci.append(locus)
    return ClinkerCluster("Query_cluster", loci)


def fasta_to_cluster(fasta_handle):
    """Convert a fasta text into a clinker.Cluster
    Args:
        fasta_handle (TextIOWrapper): handle for the fasta text
    Returns:
        a clinker.Cluster object
    """
    name = None
    start = end = 0
    sequence_length = 0
    locus_genes = []
    for line in fasta_handle:
        if line.startswith(">"):
            # if a sequence was found
            if sequence_length != 0:
                gene = ClinkerGene(label=name, start=start, end=end, strand=1)
                locus_genes.append(gene)
                # space the genes a bit
                end += FASTA_SPACE
                start = end
            name = line[1:].strip()
            sequence_length = 0
        else:
            # do not count the newline character and get in nucleotide numbers
            sequence_length += (len(line) - 1) * 3
            end += (len(line) - 1) * 3

    clinker_gene = ClinkerGene(
        label=name,
        start=start,
        end=end,
        strand=1,
        names={"accession": name}
    )
    locus_genes.append(clinker_gene)
    locus = ClinkerLocus("Locus1", locus_genes, start=0, end=end)
    return ClinkerCluster("Query_cluster", [locus])


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

    groups = {
        gene.label: ClinkerGroup(label=gene.label, genes=[gene.uid])
        for locus in clinker_query_cluster.loci
        for gene in locus.genes
    }

    for cblaster_cluster, scaffold, organism_name in cluster_hierarchies:
        clinker_cluster = cblaster_to_clinker_cluster(
            cblaster_cluster,
            scaffold.accession,
            organism_name,
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
    clinker_query_cluster = query_to_clinker_cluster(session.params["query_file"])

    # Create a Globaligner object containing mocked clusters/alignments/links
    globaligner = clusters_to_clinker_globaligner(clinker_query_cluster, cluster_hierarchies)

    if not testing:
        clinker_plot_clusters(globaligner, plot_outfile, use_file_order=True)
    if plot_outfile:
        LOG.info(f"Plot file can be found at {plot_outfile}")
    LOG.info("Done!")
