#!/usr/bin/env python3

"""
This module provides the command line interface as well as the main procedure
for running cblaster.
"""


import argparse
import logging
import sys

from cblaster import __version__, local, remote, context, helpers, database

logging.basicConfig(
    format="[%(asctime)s] %(levelname)s - %(message)s", datefmt="%H:%M:%S"
)

LOG = logging.getLogger("cblaster")
LOG.setLevel(logging.INFO)


def summarise(organisms, output=None, human=True, headers=True):
    """Generate summary of >1 Organisms, print to console or write to file.

    Parameters
    ----------
    organisms: list
    output: open file handle
    """
    summary = "\n\n\n".join(
        organism.summary(headers=headers, human=human)
        for organism in organisms
        if organism.count_hit_clusters() > 0
    )
    if summary:
        output.write(summary + "\n")
    else:
        output.write("No results found!")
    output.flush()


def count_queries(cluster, queries):
    """Count number of hits for each query in a cluster of hits."""
    counts = [0] * len(queries)
    for index, query in enumerate(queries):
        for hit in cluster:
            if query == hit.query:
                counts[index] += 1
    return counts


def value_in_args(value, args):
    """Check for value in list of argument values, removing if found.

    This function is used with arguments that take multiple values, mainly --binary.

    It checks the list for a specific value (e.g. 'hr' or 'he'). If found, the value
    is removed from the list via pop() and True is returned. If the value is not in the
    argument list, False is returned.

    This facilitates users passing in values in different orders. For example, a user
    could specify:

        --binary he path/to/filename hr

    Since this function explicitly checks for and removes 'he' and 'hr', the validation
    function can still correctly parse the filename str.
    """
    try:
        index = args.index(value)
    except ValueError:
        return False
    return args.pop(index) is not None


def validate_output_args(arguments):
    """Validate arguments for outputs.

    Checks for optional arguments 'hr' (human-readable format) and 'he'
    (show headers).
    """
    for arg in ["output", "binary"]:

        args = getattr(arguments, arg)

        # By default, want to print human-readable output with headers to stdout.
        # If binary, set flags to False and argument to None.
        if not args:
            setattr(arguments, f"{arg}_human", True if arg == "output" else False)
            setattr(arguments, f"{arg}_headers", True if arg == "output" else False)
            setattr(arguments, arg, sys.stdout if arg == "output" else None)
            continue

        setattr(arguments, f"{arg}_human", value_in_args("hr", args))
        setattr(arguments, f"{arg}_headers", value_in_args("he", args))

        if len(args) == 0:
            raise ValueError(f"No file name detected for --{arg}")

        if len(args) > 1:
            # Since we should only have a single str remaining at this point,
            # this should be sufficient validation
            raise ValueError(f"Invalid arguments provided to --{arg}")

        setattr(arguments, arg, open(args[0], "w"))


def generate_binary_table(organisms, queries, human=False, headers=True, output=None):
    """Generate a binary summary table.

    For example:

    Organism  Scaffold  Start  End    Query1  Query2  Query3  Query4
    Org 1     Scaf_1    1      20000  2       1       1       1
    Org 1     Scaf_3    3123   40302  0       1       0       1
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

    if human:
        # Calculate lengths of each column for spacing
        lengths = [max(len(row[i]) for row in rows) for i in range(columns)]
        table = "\n".join(
            "  ".join(f"{row[i]:{lengths[i]}}" for i in range(columns)) for row in rows
        )
    else:
        table = "\n".join(",".join(row) for row in rows)

    output.write(table)
    output.flush()


def makedb(genbanks, filename, indent=None):
    """Generate JSON and diamond databases."""
    db = database.DB.from_files(genbanks)

    LOG.info("Writing FASTA file with database sequences: %s", filename + ".faa")
    LOG.info("Building DIAMOND database: %s", filename + ".dmnd")
    db.makedb(filename)

    LOG.info("Building JSON database: %s", filename + ".json")
    with open(f"{filename}.json", "w") as handle:
        db.to_json(handle, indent=indent)


def cblaster(
    query_file=None,
    query_ids=None,
    mode=None,
    json=None,
    database=None,
    gap=20000,
    conserve=3,
    min_identity=30,
    min_coverage=50,
    max_evalue=0.01,
    entrez_query=None,
    output=None,
    output_human=True,
    output_headers=True,
    binary=None,
    binary_human=False,
    binary_headers=False,
    rid=None,
):
    """Run cblaster."""

    if mode == "local":
        LOG.info("Starting cblaster in local mode")
        results = local.search(
            database,
            query_file=query_file,
            query_ids=query_ids,
            min_identity=min_identity,
            min_coverage=min_coverage,
            max_evalue=max_evalue,
        )

    elif mode == "remote":
        LOG.info("Starting cblaster in remote mode")
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
    LOG.info("Fetching genomic context of hits")
    organisms = context.search(results, conserve, gap, json=json)

    if binary:
        LOG.info("Writing binary summary table to %s", binary.name)
        if query_file:
            with open(query_file) as handle:
                query_ids = list(helpers.parse_fasta(handle))

        generate_binary_table(
            organisms,
            query_ids,
            headers=binary_headers,
            human=binary_human,
            output=binary,
        )

    LOG.info("Writing summary to %s", output.name)
    summarise(organisms, output=output, human=output_human, headers=output_headers)

    return organisms


def get_arguments(args):
    """Parse arguments."""

    parser = argparse.ArgumentParser(
        "cblaster",
        description="cblaster is a tool for finding clusters of homologous"
        " proteins. Type -h/--help after either subcommand for full description of"
        " available arguments.",
        epilog="Cameron Gilchrist, 2019",
    )
    parser.add_argument(
        "--version", action="version", version="%(prog)s " + __version__
    )
    parser.add_argument(
        "-d", "--debug", help="Print debugging information", action="store_true"
    )

    subparsers = parser.add_subparsers(dest="subcommand")

    makedb = subparsers.add_parser(
        "makedb", help="Generate JSON/diamond databases from GenBank files"
    )
    makedb.add_argument(
        "genbank",
        help="Path/s to GenBank files to use when building JSON/diamond databases",
        nargs="+",
    )
    makedb.add_argument(
        "filename",
        help="Name to use when building JSON/diamond databases (with extensions"
        " .json and .dmnd, respectively)",
    )
    makedb.add_argument(
        "-i",
        "--indent",
        help="Number of spaces to use as indent in JSON database file. By default,"
        " it will be printed on a single line (indent=None) to reduce size.",
        type=int,
        default=None,
    )

    search = subparsers.add_parser(
        "search", help="Start a local/remote cblaster search"
    )

    _inputs = search.add_argument_group("Input")
    inputs = _inputs.add_mutually_exclusive_group()
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

    output = search.add_argument_group("Output")
    output.add_argument(
        "-o",
        "--output",
        nargs="+",
        help="Write results to file. Optionally, can provide 'hr' and 'he' to toggle"
        " human-readable format and headers, respectively.",
    )
    output.add_argument(
        "-b",
        "--binary",
        nargs="+",
        help="Enable output of binary table. By default, will generate a comma-"
        "delimited table; human-readable format can be specified by supplying"
        " 'hr' to this argument. Headers can be toggled by supplying 'he'."
        " e.g. cblaster --binary filename hr he ...",
    )

    searching = search.add_argument_group("Searching")
    searching.add_argument(
        "-m",
        "--mode",
        help="Mechanism through which BLAST search will be performed (def. remote)",
        choices=["local", "remote"],
        default="remote",
    )
    searching.add_argument(
        "-j",
        "--json",
        help="Path to local JSON database, created using cblaster makedb. If this"
        " argument is provided, genomic context will be fetched from this database"
        " instead of through NCBI IPG.",
    )
    searching.add_argument(
        "-db",
        "--database",
        help="Database to be searched. This should be either a path to a local"
        " DIAMOND database (if 'local' is passed to --mode) or a valid NCBI"
        " database name (def. nr)",
    )
    searching.add_argument(
        "-eq",
        "--entrez_query",
        help="An entrez search term for pre-search filtering of an NCBI database"
        " when using command line BLASTp (i.e. only used if 'remote' is passed to"
        ' --mode); e.g. "Aspergillus"[organism]',
    )
    searching.add_argument(
        "--rid",
        help="Request Identifier (RID) for a web BLAST search. This is only used"
        " if 'remote' is passed to --mode. Useful if you have previously run a web BLAST"
        " search and want to directly retrieve those results instead of running a new"
        " search.",
    )

    clusters = search.add_argument_group("Clustering")
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
        help="Minimum number of unique query sequences that must be conserved"
        " in a hit cluster for it to be reported (def. 3)",
    )
    clusters.add_argument(
        "-r",
        "--require",
        nargs="+",
        help="Names of query sequences that must be represented in a hit cluster",
    )

    filters = search.add_argument_group("Filtering")
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

    if not arguments.subcommand:
        parser.print_help()
        sys.exit()

    if arguments.subcommand == "makedb":
        return arguments

    validate_output_args(arguments)

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

    if args.subcommand == "makedb":
        makedb(args.genbank, args.filename, args.indent)
    elif args.subcommand == "search":
        cblaster(
            query_file=args.query_file,
            query_ids=args.query_ids,
            mode=args.mode,
            json=args.json,
            database=args.database,
            gap=args.gap,
            conserve=args.conserve,
            require=args.require,
            min_identity=args.min_identity,
            min_coverage=args.min_coverage,
            max_evalue=args.max_evalue,
            entrez_query=args.entrez_query,
            output=args.output,
            output_human=args.output_human,
            output_headers=args.output_headers,
            binary=args.binary,
            binary_human=args.binary_human,
            binary_headers=args.binary_headers,
            rid=args.rid,
        )


if __name__ == "__main__":
    main()
