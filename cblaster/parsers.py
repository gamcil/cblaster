"""Argument parsers."""


import argparse
import builtins

from cblaster import __version__


def add_makedb_subparser(subparsers):
    makedb = subparsers.add_parser(
        "makedb",
        help="Generate JSON/diamond databases from GenBank files"
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


def add_gui_subparser(subparsers):
    subparsers.add_parser("gui", help="Launch cblaster GUI")


def add_input_group(search):
    group = (
        search
        .add_argument_group("Input")
        .add_mutually_exclusive_group()
    )
    group.add_argument(
        "-qf",
        "--query_file",
        help="Path to FASTA file containing protein sequences to be searched",
    )
    group.add_argument(
        "-qi",
        "--query_ids",
        nargs="+",
        help="A collection of valid NCBI sequence identifiers to be searched",
    )


def add_output_arguments(group):
    group.add_argument(
        "-o",
        "--output",
        help="Write results to file",
    )
    group.add_argument(
        "-ohh",
        "--output_hide_headers",
        action="store_true",
        help="Hide headers when printing result output."
    )
    group.add_argument(
        "-ode",
        "--output_delimiter",
        help="Delimiter character to use when printing result output.",
        default=None,
    )
    group.add_argument(
        "-odc",
        "--output_decimals",
        type=int,
        help="Total decimal places to use when printing score values",
        default=4,
    )


def add_binary_arguments(group):
    group.add_argument(
        "-b",
        "--binary",
        help="Generate a binary table.",
    )
    group.add_argument(
        "-bhh",
        "--binary_hide_headers",
        action="store_true",
        help="Hide headers in the binary table.",
    )
    group.add_argument(
        "-bde",
        "--binary_delimiter",
        help="Delimiter used in binary table (def. none = human readable).",
        default=None,
    )
    group.add_argument(
        "-bkey",
        "--binary_key",
        help="Key function used when generating binary table cell values.",
        default="len",
        choices=["len", "max", "sum"],
    )
    group.add_argument(
        "-bat",
        "--binary_attr",
        help="Hit attribute used when generating binary table cell values.",
        default="identity",
        choices=["identity", "coverage", "bitscore", "evalue"],
    )
    group.add_argument(
        "-bdc",
        "--binary_decimals",
        help="Total decimal places to use when printing score values",
        default=4,
    )


def add_output_group(search):
    group = search.add_argument_group("Output")

    add_output_arguments(group)
    add_binary_arguments(group)

    group.add_argument(
        "-p",
        "--plot",
        nargs="?",
        const=True,
        default=False,
        help="Generate a cblaster plot. If this argument is specified with no"
        " file name, the plot will be served using Python's HTTP server. If a"
        " file name is specified, a static HTML file will be generated at that"
        " path."
    )


def add_searching_group(search):
    group = search.add_argument_group("Searching")
    group.add_argument(
        "-m",
        "--mode",
        help="cblaster search mode",
        choices=["local", "remote"],
        default="remote",
    )
    group.add_argument(
        "-db",
        "--database",
        default="nr",
        help="Database to be searched. This should be either a path to a local"
        " DIAMOND database (if 'local' is passed to --mode) or a valid NCBI"
        " database name (def. nr)",
    )
    group.add_argument(
        "-jdb",
        "--json_db",
        help="Path to local JSON database created using cblaster makedb. If this"
        " argument is provided, genomic context will be fetched from this database"
        " instead of through NCBI IPG.",
    )
    group.add_argument(
        "-eq",
        "--entrez_query",
        help="An NCBI Entrez search term for pre-search filtering of an NCBI database"
        " when using command line BLASTp (i.e. only used if 'remote' is passed to"
        ' --mode); e.g. "Aspergillus"[organism]',
    )
    group.add_argument(
        "--rid",
        help="Request Identifier (RID) for a web BLAST search. This is only used"
        " if 'remote' is passed to --mode. Useful if you have previously run a web BLAST"
        " search and want to directly retrieve those results instead of running a new"
        " search.",
    )
    group.add_argument(
        "-s",
        "--session_file",
        help="Load session from JSON. If the specified file does not exist, "
        "the results of the new search will be saved to this file.",
    )
    group.add_argument(
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


def add_clustering_group(search):
    group = search.add_argument_group("Clustering")
    group.add_argument(
        "-g",
        "--gap",
        type=int,
        default=20000,
        help="Maximum allowed intergenic distance (bp) between conserved hits to"
        " be considered in the same block (def. 20000)",
    )
    group.add_argument(
        "-u",
        "--unique",
        type=int,
        default=3,
        help="Minimum number of unique query sequences that must be conserved"
        " in a hit cluster (def. 3)",
    )
    group.add_argument(
        "-mh",
        "--min_hits",
        type=int,
        default=3,
        help="Minimum number of hits in a cluster (def. 3)",
    )
    group.add_argument(
        "-r",
        "--require",
        nargs="+",
        help="Names of query sequences that must be represented in a hit cluster",
    )


def add_filtering_group(search):
    group = search.add_argument_group("Filtering")
    group.add_argument(
        "-me",
        "--max_evalue",
        type=float,
        default=0.01,
        help="Maximum e-value for a BLAST hit to be saved (def. 0.01)",
    )
    group.add_argument(
        "-mi",
        "--min_identity",
        type=float,
        default=30,
        help="Minimum percent identity for a BLAST hit to be saved (def. 30)",
    )
    group.add_argument(
        "-mc",
        "--min_coverage",
        type=float,
        default=50,
        help="Minimum percent query coverage for a BLAST hit to be saved (def. 50)",
    )


def add_search_subparser(subparsers):
    search = subparsers.add_parser("search", help="Start a local/remote cblaster search")
    add_input_group(search)
    add_output_group(search)
    add_searching_group(search)
    add_clustering_group(search)
    add_filtering_group(search)


def add_gne_subparser(subparsers):
    gne = subparsers.add_parser("gne", help="Perform gene neighbourhood estimation")
    gne.add_argument("session", help="cblaster session file")
    gne.add_argument("--output", help="Path to output file")
    gne.add_argument("--delimiter", default=",", help="Character used as delimiter")
    gne.add_argument("--max_gap", type=int, default=1000000, help="Maximum intergenic distance")
    gne.add_argument("--samples", type=int, default=1000, help="Total samples taken from max_gap")
    gne.add_argument(
        "--scale",
        choices=["linear", "log"],
        default="log",
        help="Draw sampling values from a linear or log scale"
    )


def get_parser():
    parser = argparse.ArgumentParser(
        "cblaster",
        description="cblaster is a tool for finding clusters of homologous"
        " proteins. Type -h/--help after either subcommand for full description of"
        " available arguments.",
        epilog="Cameron Gilchrist, 2020",
    )
    parser.add_argument("--version", action="version", version="%(prog)s " + __version__)
    parser.add_argument("-d", "--debug", help="Print debugging information", action="store_true")
    parser.add_argument(
        "-i",
        "--indent",
        help="Total spaces to use as indent in JSON file (def. None)",
        type=int,
        default=None,
    )
    subparsers = parser.add_subparsers(dest="subcommand")
    add_gui_subparser(subparsers)
    add_makedb_subparser(subparsers)
    add_search_subparser(subparsers)
    add_gne_subparser(subparsers)
    return parser


def parse_args(args):
    parser = get_parser()
    arguments = parser.parse_args(args)

    if not arguments.subcommand:
        parser.print_help()
        raise SystemExit

    if arguments.subcommand in ("gui", "makedb", "gne"):
        return arguments

    if arguments.mode == "remote":
        valid_dbs = ("nr", "refseq_protein", "swissprot", "pdbaa")
        if arguments.database not in valid_dbs:
            parser.error(f"Valid databases are: {', '.join(valid_dbs)}")
    else:
        for arg in ["entrez_query", "rid"]:
            if getattr(arguments, arg):
                parser.error(f"--{arg} can only be used when --mode is 'remote'")

    # Convert key to its corresponding builtin function
    arguments.binary_key = getattr(builtins, arguments.binary_key)

    if arguments.recompute and not arguments.session_file:
        parser.error("--recompute requires --session_file")

    return arguments
