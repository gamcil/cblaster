#!/usr/bin/env python3

"""
This module provides the command line interface as well as the main procedure
for running cblaster.
"""


import argparse
import logging
import sys
from pathlib import Path

from cblaster import __version__, context, database, helpers, local, plot, remote
from cblaster.classes import Session

logging.basicConfig(
    format="[%(asctime)s] %(levelname)s - %(message)s", datefmt="%H:%M:%S"
)

LOG = logging.getLogger("cblaster")
LOG.setLevel(logging.INFO)


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


def makedb(genbanks, filename, indent=None):
    """Generate JSON and diamond databases."""
    db = database.DB.from_files(genbanks)

    LOG.info("Writing FASTA file with database sequences: %s", filename + ".faa")
    LOG.info("Building DIAMOND database: %s", filename + ".dmnd")
    db.makedb(filename)

    LOG.info("Building JSON database: %s", filename + ".json")
    with open(f"{filename}.json", "w") as handle:
        db.to_json(handle, indent=indent)


def filter_session(
    session, min_identity, min_coverage, max_evalue, gap, unique, min_hits, require
):
    """Filter a previous session with new thresholds."""

    def hit_meets_thresholds(hit):
        return (
            hit.identity > min_identity
            and hit.coverage > min_coverage
            and hit.evalue < max_evalue
        )

    def filter_scaffold(scaffold):
        scaffold.hits = [hit for hit in scaffold.hits if hit_meets_thresholds(hit)]
        scaffold.clusters = list(
            context.find_clusters(
                scaffold.hits,
                unique=unique,
                min_hits=min_hits,
                gap=gap,
                require=require,
            )
        )

    for organism in session.organisms:
        for scaffold in organism.scaffolds.values():
            filter_scaffold(scaffold)


def cblaster(
    query_file=None,
    query_ids=None,
    mode=None,
    json_db=None,
    database=None,
    gap=20000,
    unique=3,
    min_hits=3,
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
    require=None,
    session_file=None,
    indent=None,
    figure=False,
    figure_dpi=300,
    recompute=False,
):
    """Run cblaster."""

    if session_file and Path(session_file).exists():
        LOG.info("Loading %s", session_file)
        with open(session_file) as fp:
            session = Session.from_json(fp)

        if recompute:
            LOG.info("Filtering session with new thresholds")
            filter_session(
                session,
                min_identity,
                min_coverage,
                max_evalue,
                gap,
                unique,
                min_hits,
                require,
            )

            if recompute is not True:
                LOG.info("Writing recomputed session to %s", recompute)
                with open(recompute, "w") as fp:
                    session.to_json(fp, indent=indent)
    else:
        session = Session(
            query_ids if query_ids else [],
            params={
                "mode": mode,
                "database": database,
                "min_identity": min_identity,
                "min_coverage": min_coverage,
                "max_evalue": max_evalue,
            },
        )

        if query_file:
            with open(query_file) as fp:
                sequences = helpers.parse_fasta(fp)
            session.queries = list(sequences)
            session.params["query_file"] = query_file

        if json_db:
            session.params["json_db"] = json_db

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

            if entrez_query:
                session.params["entrez_query"] = entrez_query

            rid, results = remote.search(
                query_file=query_file,
                query_ids=query_ids,
                rid=rid,
                database=database,
                min_identity=min_identity,
                min_coverage=min_coverage,
                max_evalue=max_evalue,
                entrez_query=entrez_query,
            )

            session.params["rid"] = rid

        LOG.info("Found %i hits meeting score thresholds", len(results))
        LOG.info("Fetching genomic context of hits")
        session.organisms = context.search(
            results,
            unique=unique,
            min_hits=min_hits,
            gap=gap,
            require=require,
            json_db=json_db,
        )

        if session_file:
            LOG.info("Writing current search session to %s", session_file)
            with open(session_file, "w") as fp:
                session.to_json(fp, indent=indent)

    if binary:
        LOG.info("Writing binary summary table to %s", binary.name)
        session.format("binary", binary, human=binary_human, headers=binary_headers)

    LOG.info("Writing summary to %s", output.name)
    session.format("summary", output, human=output_human, headers=output_headers)

    if figure:
        if figure is True:
            LOG.info("Generating cblaster plot...")
            plot.plot(session)
        else:
            LOG.info("Writing figure to %s", figure)
            plot.plot(session, figure=figure, dpi=figure_dpi)

    return session


def get_arguments(args):
    """Parse arguments."""

    parser = argparse.ArgumentParser(
        "cblaster",
        description="cblaster is a tool for finding clusters of homologous"
        " proteins. Type -h/--help after either subcommand for full description of"
        " available arguments.",
        epilog="Cameron Gilchrist, 2020",
    )
    parser.add_argument(
        "--version", action="version", version="%(prog)s " + __version__
    )
    parser.add_argument(
        "-d", "--debug", help="Print debugging information", action="store_true"
    )
    parser.add_argument(
        "-i",
        "--indent",
        help="Total spaces to use as indent in JSON file (def. None)",
        type=int,
        default=None,
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
    output.add_argument(
        "-f",
        "--figure",
        nargs="?",
        const=True,
        default=False,
        help="Save a cblaster plot to file. File type is determined from the file"
        " extension (e.g. --figure x.svg will save as SVG). If no file is specified,"
        " the interactive matplotlib viewer will be opened.",
    )
    output.add_argument(
        "-dpi",
        "--figure_dpi",
        help="DPI to use when saving figure (def. 300)",
        default=300,
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
        "-jdb",
        "--json_db",
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
    searching.add_argument(
        "-s",
        "--session_file",
        help="Load session from JSON. If the specified file does not exist, "
        "the results of the new search will be saved to this file.",
    )
    searching.add_argument(
        "-rcp",
        "--recompute",
        nargs="?",
        const=True,
        default=False,
        help="Recompute previous search session using new thresholds. The filtered"
        " session will be written to the file specified by this argument. If this"
        " argument is specified with no value, the session will be filtered but"
        " not saved (e.g. for plotting purposes).",
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
        "-u",
        "--unique",
        type=int,
        default=3,
        help="Minimum number of unique query sequences that must be conserved"
        " in a hit cluster (def. 3)",
    )
    clusters.add_argument(
        "-mh",
        "--min_hits",
        type=int,
        default=3,
        help="Minimum number of hits in a cluster (def. 3)",
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
            parser.error(f"Valid databases are: {', '.join(valid_dbs)}")
    else:
        for arg in ["entrez_query", "rid"]:
            if getattr(arguments, arg):
                parser.error(f"--{arg} can only be used when --mode is 'remote'")

    if arguments.recompute and not arguments.session_file:
        parser.error("-rcp/--recompute requires -s/--session_file")

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
            json_db=args.json_db,
            database=args.database,
            gap=args.gap,
            unique=args.unique,
            min_hits=args.min_hits,
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
            session_file=args.session_file,
            indent=args.indent,
            recompute=args.recompute,
            figure_dpi=args.figure_dpi,
            figure=args.figure,
        )


if __name__ == "__main__":
    main()
