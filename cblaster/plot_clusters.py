
import tempfile
from pathlib import Path
import subprocess
import os
import shutil
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


def find_genbank_files(files):
    genbank_files = []
    for path in files:
        path_obj = Path(path)
        if path_obj.is_dir():
            genbank_files.extend(
                [str(po.resolve()) for po in path_obj.iterdir() if po.suffix in (".gbk", ".gb", ".genbank", ".gbff")])
        elif path_obj.suffix in (".gbk", ".gb", ".genbank", ".gbff"):
            genbank_files.append(str(path_obj.resolve()))
    return genbank_files


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


def clusters_to_clinker_allignments(query_cluster, both_clusters):
    allignments = []
    for cblaster_cluster, clinker_cluster in both_clusters:
        allignment = ClinkerAlignment(query=query_cluster, target=clinker_cluster)
        for subject in cblaster_cluster.subjects:
            best_hit = max(subject.hits, key=lambda x: x.bitscore)
            query_gene = _gene_from_clinker_cluster(query_cluster, best_hit.query)
            subject_gene = _gene_from_clinker_cluster(clinker_cluster, best_hit.subject)
            allignment.add_link(query_gene, subject_gene, best_hit.identity, 0)
        allignments.append(allignment)
    return allignments


def _gene_from_clinker_cluster(cluster, gene_label):
    """Here because of a bug in clinker where a name attribute is called when that should be label"""
    for locus in cluster.loci:
        for gene in locus.genes:
            if gene.label == gene_label:
                return gene


def allignments_to_clinker_global_alligner(allignments):
    global_aligner = ClinkerGlobalaligner()
    for allignment in allignments:
        global_aligner.add_alignment(allignment)
    return global_aligner


def plot_clusters(
        session=None,
        files=None,
        cluster_numbers=None,
        score_threshold=None,
        organisms=None,
        scaffolds=None,
        allign_clusters=False,
        identity=0.3,
        plot_outfile=None,
        allignment_out=None,
        cluster_out=None,
        prefix="",
):
    with open(session, "r") as f:
        session = Session.from_json(f.read())
    cluster_hierarchies = extract_cluster_hierarchies(session, cluster_numbers, score_threshold, organisms, scaffolds)
    both_clusters = []
    for cluster, scaffold_acs, org_name in cluster_hierarchies:
        both_clusters.append((cluster, cluster.to_clinker_cluster()))
    clinker_query_cluster = query_to_clinker_cluster(session.params["query_file"])

    allignments = clusters_to_clinker_allignments(clinker_query_cluster, both_clusters)
    global_aligner = allignments_to_clinker_global_alligner(allignments)

    clinker_plot_clusters(global_aligner, plot_outfile, use_file_order=True)
    LOG.info(f"Plot file can be found at {plot_outfile}")
    LOG.info("Done!")



# def plot_clusters(
#     session=None,
#     files=None,
#     cluster_numbers=None,
#     score_threshold=None,
#     organisms=None,
#     scaffolds=None,
#     allign_clusters=False,
#     identity=0.3,
#     plot_outfile=None,
#     allignment_out=None,
#     cluster_out=None,
#     prefix="",
# ):
#     # if no genbank files are provided make sure to create them
#     remove_temp = False
#     if not files:
#         if not cluster_out:
#             cluster_out = tempfile.mkdtemp()
#             remove_temp = True
#         extract_clusters(
#             session,
#             cluster_out,
#             prefix=prefix,
#             cluster_numbers=cluster_numbers,
#             score_threshold=score_threshold,
#             organisms=organisms,
#             scaffolds=scaffolds,
#         )
#         files = [cluster_out]
#
#     # get only the genbank files present in directories and files
#     genbank_files = find_genbank_files(files)
#     # makes sure to add the query as well for comparisson
#     with open(session, "r") as f:
#         session = Session.from_json(f.read())
#         genbank_files.append(session.params["query_file"])
#
#     try:
#         run_clinker(genbank_files, allign_clusters, identity, plot_outfile, allignment_out)
#     except subprocess.CalledProcessError:
#         # make sure to remove the temp dir even when clinker crashes
#         # if remove_temp:
#         #     shutil.rmtree(cluster_out)
#         raise SystemExit
#     # make sure to remove the temp dir
#     if remove_temp:
#         shutil.rmtree(cluster_out)
#     LOG.info(f"Plot file can be found at {plot_outfile}")
#     LOG.info("Done!")
