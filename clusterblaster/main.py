#!/usr/bin/env python3

"""
This module provides the command line interface as well as the main procedure
for running clusterblaster.
"""


import argparse
import logging
import sys

from clusterblaster import __version__, local, remote, context, helpers

logging.basicConfig(
    format="[%(asctime)s] %(levelname)s - %(message)s", datefmt="%H:%M:%S"
)

LOG = logging.getLogger("clusterblaster")
LOG.setLevel(logging.INFO)


def summarise(organisms, output=None):
    """Generate summary of >1 Organisms, print to console or write to file.

    Parameters
    ----------
    organisms: list
    output: open file handle
    """
    summary = "\n\n\n".join(
        organism.summary()
        for organism in organisms
        if organism.count_hit_clusters() > 0
    )
    if output:
        output.write(summary)
    else:
        print(summary, flush=True)


def count_queries(cluster, queries):
    """Count number of hits for each query in a cluster of hits."""
    counts = [0] * len(queries)
    for index, query in enumerate(queries):
        for hit in cluster:
            if query == hit.query:
                counts[index] += 1
    return counts


def generate_binary_table(organisms, queries, headers=True, output=None):
    """Generate a binary summary table.

    Format:
    Organism  Scaffold  Query1  Query2  Query3  Query4
    Org 1     Scaf_1    2       1       1       1
    Org 1     Scaf_3    0       1       0       1
    Org 2     Scaf_1    1       0       1       1
    """
    columns = len(queries) + 4

    rows = [
        [
            organism.full_name,
            accession,
            str(cluster[0].start),
            str(cluster[-1].end),
            *(str(count) for count in count_queries(cluster, queries)),
        ]
        for organism in organisms
        for accession, scaffold in organism.scaffolds.items()
        for cluster in scaffold.clusters
    ]

    if headers:
        rows.insert(0, ["Organism", "Scaffold", "Start", "End", *queries])

    lengths = [max(len(row[i]) for row in rows) for i in range(columns)]

    table = "\n".join(
        "  ".join(f"{row[i]:{lengths[i]}}" for i in range(columns)) for row in rows
    )

    if output:
        output.write(table)
    else:
        print(table, flush=True)


def clusterblaster(
    query_file=None,
    query_ids=None,
    mode=None,
    database=None,
    gap=20000,
    conserve=3,
    min_identity=30,
    min_coverage=50,
    max_evalue=0.01,
    entrez_query=None,
    output=None,
    binary=None,
    rid=None,
):
    """Run clusterblaster."""

    if mode == "local":
        LOG.info("Starting clusterblaster in local mode")
        results = local.search(
            database,
            query_file=query_file,
            query_ids=query_ids,
            min_identity=min_identity,
            min_coverage=min_coverage,
            max_evalue=max_evalue,
            entrez_query=entrez_query,
        )

    elif mode == "remote":
        LOG.info("Starting clusterblaster in remote mode")
        results = remote.search(
            query_file=query_file,
            query_ids=query_ids,
            rid=rid,
            database=database,
            min_identity=min_identity,
            min_coverage=min_coverage,
            max_evalue=max_evalue,
            entrez_query=entrez_query,
        )

    LOG.info("Found %i hits meeting score thresholds", len(results))
    LOG.info("Fetching genomic context of hits from NCBI")
    organisms = context.search(results, conserve, gap)

    if binary:
        LOG.info("Writing binary summary table to %s", binary.name)
        if query_file:
            with open(query_file) as handle:
                query_ids = list(helpers.parse_fasta(handle))

        generate_binary_table(organisms, query_ids, output=binary)

    if output:
        LOG.info("Writing summary to %s", output.name)
        summarise(organisms, output=output)
    else:
        summarise(organisms)

    return organisms


def get_arguments(args):
    """Parse arguments."""

    parser = argparse.ArgumentParser(
        "clusterblaster",
        description="clusterblaster is a tool for finding clusters of homologous"
        " proteins.",
        epilog="Cameron Gilchrist, 2019",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=argparse.FileType("w"),
        help="Save output to this file path",
    )
    parser.add_argument(
        "-b",
        "--binary",
        type=argparse.FileType("w"),
        help="Save binary format table to this file path",
    )
    parser.add_argument(
        "--version", action="version", version="%(prog)s " + __version__
    )
    parser.add_argument(
        "-d", "--debug", help="Print debugging information", action="store_true"
    )

    _inputs = parser.add_argument_group("Input")
    inputs = _inputs.add_mutually_exclusive_group(required=True)
    inputs.add_argument(
        "-qf",
        "--query_file",
        help="Path to FASTA file containing protein sequences to be searched",
    )
    inputs.add_argument(
        "-qi",
        "--query_ids",
        nargs="+",
        help="A collection of valid NCBI sequence identifiers to be searched",
    )

    search = parser.add_argument_group("Searching")
    search.add_argument(
        "-m",
        "--mode",
        help="Mechanism through which BLAST search will be performed (def. remote)",
        choices=["local", "remote"],
        default="remote",
    )
    search.add_argument(
        "-db",
        "--database",
        help="Database to be searched. This should be either a path to a local"
        " DIAMOND database (if 'local' is passed to --mode) or a valid NCBI"
        " database name (def. nr)",
    )
    search.add_argument(
        "-eq",
        "--entrez_query",
        help="An entrez search term for pre-search filtering of an NCBI database"
        " when using command line BLASTp (i.e. only used if 'remote' is passed to"
        ' --mode); e.g. "Aspergillus"[organism]',
    )
    search.add_argument(
        "--rid",
        help="Request Identifier (RID) for a web BLAST search. This is only used"
        " if 'remote' is passed to --mode. Useful if you have previously run a web BLAST"
        " search and want to directly retrieve those results instead of running a new"
        " search.",
    )

    clusters = parser.add_argument_group("Clustering")
    clusters.add_argument(
        "-g",
        "--gap",
        type=int,
        default=20000,
        help="Maximum allowed intergenic distance (bp) between conserved hits to"
        " be considered in the same block (def. 20000)",
    )
    clusters.add_argument(
        "-c",
        "--conserve",
        type=int,
        default=3,
        help="Minimum number of query sequences that must be conserved in a single"
        " block for it to be reported (def. 3)",
    )

    filters = parser.add_argument_group("Filtering")
    filters.add_argument(
        "-me",
        "--max_evalue",
        type=float,
        default=0.01,
        help="Maximum e-value for a BLAST hit to be saved (def. 0.01)",
    )
    filters.add_argument(
        "-mi",
        "--min_identity",
        type=float,
        default=30,
        help="Minimum percent identity for a BLAST hit to be saved (def. 30)",
    )
    filters.add_argument(
        "-mc",
        "--min_coverage",
        type=float,
        default=50,
        help="Minimum percent query coverage for a BLAST hit to be saved (def. 50)",
    )

    arguments = parser.parse_args(args)

    if arguments.mode == "remote":
        if not arguments.database:
            # Default to non-redundant database if --database is not given
            arguments.database = "nr"

        valid_dbs = ("nr", "refseq_protein", "swissprot", "pdbaa")
        if arguments.database not in valid_dbs:
            raise ValueError(f"Valid databases are: {', '.join(valid_dbs)}")
    else:
        for arg in ["entrez_query", "rid"]:
            if getattr(arguments, arg):
                raise ValueError(f"--{arg} can only be used when --mode is 'remote'")

    return arguments


def main():
    args = get_arguments(sys.argv[1:])

    if args.debug:
        LOG.setLevel(logging.DEBUG)

    clusterblaster(
        query_file=args.query_file,
        query_ids=args.query_ids,
        mode=args.mode,
        database=args.database,
        gap=args.gap,
        conserve=args.conserve,
        min_identity=args.min_identity,
        min_coverage=args.min_coverage,
        max_evalue=args.max_evalue,
        entrez_query=args.entrez_query,
        output=args.output,
        binary=args.binary,
        rid=args.rid,
    )


if __name__ == "__main__":
    main()
