"""Argument parsers."""


import argparse
import builtins
from pathlib import Path
import os
from multiprocessing import cpu_count

from cblaster import __version__


NCBI_DATABASES = ("nr", "refseq_protein", "swissprot", "pdbaa")


def full_path(file_path, *acces_modes, dir=False):
    """Test if a file path or directory exists and has the correct permissions and create a full path

    For reading acces the file has to be pressent and there has to be read acces. For writing acces the directory with
    the file has to be present and there has to be write acces in that directory.

    Args:
        file_path (str): relative or absoluete path to a file
        acces_modes (List): a list of integers of acces modes for which at least one should be allowed
        dir (bool): if the path is to a directory or not
    Returns:
        A string that is the full path to the provided file_path
    Raises:
        argparse.ArgumentTypeError when the provided path does not exist or the file does not have the correct
        permissions to be accessed
    """
    full_file_path = Path(file_path).absolute().resolve()
    failed_path = failed_acces = False
    for acces_mode in acces_modes:
        if not dir and full_file_path.is_file() or (acces_mode == os.W_OK and full_file_path.parent.is_dir()):
            if os.access(full_file_path, acces_mode) or \
                    (acces_mode == os.W_OK and os.access(full_file_path.parent, acces_mode)):
                return str(full_file_path)
            else:
                failed_acces = True
        elif dir and full_file_path.is_dir():
            if os.access(full_file_path, acces_mode):
                return str(full_file_path)
        else:
            failed_path = True
    if failed_path:
        raise argparse.ArgumentTypeError(f"Invalid path: '{file_path}'.")
    elif failed_acces:
        raise argparse.ArgumentTypeError(f"Invalid acces for path: {file_path}.")


def full_database_path(database, *acces_modes):
    """Make sure the database path is also correct, but do not check when providing one of the NCBI databases

    Args:
        database (str): a string that is the path to the database creation files or a NCBI database identifier
        acces_modes (List): a list of integers of acces modes for which at least one should be allowed
    Returns:
        a string that is the full path to the database file or a NCBI database identifier
    """
    if database not in NCBI_DATABASES:
        try:
            return full_path(database, *acces_modes)
        except argparse.ArgumentTypeError as e:
            raise type(e)(str(e) + f" Or use one of the following databases {', '.join(NCBI_DATABASES)}"
                                   f" when running in remote mode")
    return database


def max_cpus(value):
    """
    Ensure that the cpu's do not go above the available amount. Setting to high cpu's will crash database creation
    badly

    Args:
        value(int): number of cpu's as provided by the user
    Returns:
        value as an integer with 1 <= value <= multiprocessing.cpu_count()
    """
    try:
        value = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("Invalid 'int' value: f")
    value = max(1, min(value, cpu_count()))
    return value


def add_makedb_subparser(subparsers):
    description = "Generate local databases from genome files"
    makedb = subparsers.add_parser("makedb", help=description, description=description)
    makedb.add_argument(
        "paths",
        type=lambda x: full_path(x, os.R_OK),
        help="Path/s to genome files to use when building local databases",
        nargs="+",
    )
    makedb.add_argument(
        "-n",
        "--name",
        required=True,
        help="Name to use when building sqlite3/diamond databases (with extensions"
             " .sqlite3 and .dmnd, respectively)",
    )
    makedb.add_argument(
        "-cp",
        "--cpus",
        type=max_cpus,
        help="Number of CPUs to use when parsing genome files. By default, all"
             " available cores will be used.",
    )
    makedb.add_argument(
        "-b",
        "--batch",
        type=int,
        help="Number of genome files to parse before saving them in the local"
             " database. Useful when encountering memory issues with large/many"
             " files. By default, all genome files will be parsed at once."
    )
    makedb.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="Overwrite pre-existing files, if any"
    )


def add_gui_subparser(subparsers):
    subparsers.add_parser("gui", help="Launch cblaster GUI")


def add_input_group(search):
    group = search.add_argument_group("Input")
    group.add_argument(
        "-qf",
        "--query_file",
        type=lambda x: full_path(x, os.R_OK),
        help="Path to FASTA file containing protein sequences to be searched",
    )
    group.add_argument(
        "-qi",
        "--query_ids",
        nargs="+",
        help="A collection of valid NCBI sequence identifiers to be searched",
    )
    group.add_argument(
        "-qp",
        "--query_profiles",
        nargs="+",
        help="A collection of valid Pfam profile identifiers to be searched "
    )


def add_output_arguments(group):
    group.add_argument(
        "-o",
        "--output",
        type=lambda x: full_path(x, os.W_OK),
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
    group.add_argument(
        "-osc",
        "--sort_clusters",
        action="store_true",
        help="Sorts the clusters of the final output on score. This means that clusters of the same organism are not"
             " neccesairily close together in the output."
    )


def add_binary_arguments(group):
    group.add_argument(
        "-b",
        "--binary",
        type=lambda x: full_path(x, os.W_OK),
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
        type=lambda x: full_path(x, os.W_OK),
        help="Generate a cblaster plot. If this argument is specified with no"
             " file name, the plot will be served using Python's HTTP server. If a"
             " file name is specified, a static HTML file will be generated at that"
             " path."
    )
    group.add_argument(
        "-mpc",
        "--max_plot_clusters",
        default=50,
        type=int,
        help="The maximum amount of clusters included in the plot when sorting clusters on"
             " score, meaning -osc has to be used for this argument to take effect. (def 50)"
    )
    group.add_argument(
        "--blast_file",
        type=lambda x: full_path(x, os.W_OK),
        help="Save BLAST/DIAMOND hit table to file"
    )
    group.add_argument(
        "--ipg_file",
        type=lambda x: full_path(x, os.W_OK),
        help="Save IPG table to file (only if --mode remote)"
    )


def add_searching_group(search):
    group = search.add_argument_group("Searching")
    group.add_argument(
        "-m",
        "--mode",
        help="cblaster search mode",
        choices=["local", "remote", "hmm", "combi_local", "combi_remote"],
        default="remote",
    )
    group.add_argument(
        "-db",
        "--database",
        default=["nr"],
        nargs="+",
        type=lambda x: full_database_path(x, os.R_OK),
        help="Database to be searched. Remote search mode: NCBI database name (def."
        " 'nr'); local search mode: path to DIAMOND database; HMM search mode:"
        " path to FASTA file. In local/hmm/combined modes, must have cblaster"
        " database in same location with same name and .sqlite3 extension."
    )
    group.add_argument(
        "-cp",
        "--cpus",
        type=int,
        help="Number of CPUs to use in local search. By default, all"
             " available cores will be used.",
    )
    group.add_argument(
        "-pfam",
        "--database_pfam",
        help="Path to folder containing Pfam database files (Pfam-A.hmm.gz and"
        " Pfam-A.dat.gz). If not found, cblaster will download the latest Pfam"
        " release to this folder. This option is required when running HMM or"
        " combi search modes.",
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
        nargs="*",
        type=lambda x: full_path(x, os.R_OK, os.W_OK),
        help="Load session from JSON. If the specified file does not exist, "
             "the results of the new search will be saved to this file.",
    )
    group.add_argument(
        "-rcp",
        "--recompute",
        nargs="?",
        const=True,
        default=False,
        type=lambda x: full_path(x, os.W_OK),
        help="Recompute previous search session using new thresholds. The filtered"
             " session will be written to the file specified by this argument. If this"
             " argument is specified with no value, the session will be filtered but"
             " not saved (e.g. for plotting purposes).",
    )
    group.add_argument(
        "-hs",
        "--hitlist_size",
        type=int,
        default=5000,
        help="Maximum total hits to save from a remote BLAST search (def. 5000). Setting"
             " this value too low may result in missed hits/clusters."
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
    group.add_argument(
        "-per",
        "--percentage",
        type=int,
        default=50,
        help="Filter on %% of query genes needed to be present in cluster",
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


def add_intermediate_genes_group(search):
    group = search.add_argument_group("Intermediate genes")

    group.add_argument(
        "-ig",
        "--intermediate_genes",
        action="store_true",
        help="Show genes that in or near clusters but not part of the cluster. "
             "This takes some extra computation time."
    )
    group.add_argument(
        "-md",
        "--max_distance",
        type=int,
        default=5000,
        help="The maximum distance between the start/end of a cluster and an intermediate gene (def. 5000)"
    )
    group.add_argument(
        "-mic",
        "--maximum_clusters",
        type=int,
        default=100,
        help="The maximum amount of clusters will get intermediate genes assigned. Ordered on score (def. 100)"
    )


def add_search_subparser(subparsers):
    search = subparsers.add_parser(
        "search",
        help="Start a local/remote cblaster search",
        description="Remote/local cblaster searches.",
        epilog="Example usage\n-------------\n"
               "Run a remote cblaster search, save the session and generate a plot:\n"
               "  $ cblaster search -qf query.fa -s session.json -p\n\n"
               "Recompute a search session with new parameters:\n"
               "  $ cblaster search -s session.json -rcp new.json -u 4 -g 40000\n\n"
               "Merge multiple search sessions:\n"
               "  $ cblaster search -s one.json two.json three.json -rcp merged.json\n\n"
               "Perform a local search:\n"
               "  $ cblaster makedb $(ls folder/*.gbk) mydb\n"
               "  $ cblaster search -qf query.fa -db mydb.dmnd -jdb mydb.json\n\n"
               "Save plot as a static HTML file:\n"
               "  $ cblaster search -s session.json -p gne.html\n\n"
               "Kitchen sink example:\n"
               "  $ cblaster search --query_file query.fa \\ \n"
               "      --session_file session.json \\ \n"
               "      --plot my_plot.html \\ \n"
               "      --output summary.csv --output_decimals 2 \\ \n"
               "      --binary abspres.csv --binary_delimiter \",\" \\ \n"
               "      --entrez_query \"Aspergillus\"[orgn] \\ \n"
               "      --max_evalue 0.05 --min_identity 50 --min_coverage 70 \\ \n"
               "      --gap 50000 --unique 2 --min_hits 3 --require Gene1 Gene2\n\n"
               "Cameron Gilchrist, 2020",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    add_input_group(search)
    add_output_group(search)
    add_searching_group(search)
    add_clustering_group(search)
    add_filtering_group(search)
    add_intermediate_genes_group(search)


def add_gne_output_group(parser):
    group = parser.add_argument_group("Output")
    group.add_argument(
        "-o",
        "--output",
        type=lambda x: full_path(x, os.W_OK),
        help="Write results to file",
    )
    group.add_argument(
        "-hh",
        "--hide_headers",
        action="store_true",
        help="Hide headers when printing result output."
    )
    group.add_argument(
        "-d",
        "--delimiter",
        help="Delimiter character to use when printing result output.",
        default=None,
    )
    group.add_argument(
        "-e",
        "--decimals",
        type=int,
        help="Total decimal places to use when printing score values",
        default=4,
    )
    group.add_argument(
        "-p",
        "--plot",
        nargs='?',
        const=True,
        default=None,
        type=lambda x: full_path(x, os.W_OK),
        help="Specify this argument without value to dynamically serve te plot. If a file location is provided"
             " the plot will be saved there."
    )


def add_gne_params_group(parser):
    group = parser.add_argument_group("Parameters")
    group.add_argument(
        "--max_gap",
        type=int,
        default=100000,
        help="Maximum intergenic distance (def. 100000)"
    )
    group.add_argument(
        "--samples",
        type=int,
        default=100,
        help="Total samples taken from max_gap (def. 100)"
    )
    group.add_argument(
        "--scale",
        choices=["linear", "log"],
        default="linear",
        help="Draw sampling values from a linear or log scale (def. linear)"
    )


def add_gne_subparser(subparsers):
    gne = subparsers.add_parser(
        "gne",
        help="Perform gene neighbourhood estimation",
        description="Gene neighbourhood estimation.\n"
                    "Repeatedly recomputes homologue clusters with different --gap values.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Example usage\n-------------\n"
               "Maximum gap value 200Kbp, with 200 evenly distributed gap values:\n"
               "  $ cblaster gne session.json --max_gap 200000 --samples 200 --scale linear\n\n"
               "Draw gap values from a log scale (gaps increase as values increase):\n"
               "  $ cblaster gne session.json --scale log\n\n"
               "Save delimited tabular output:\n"
               "  $ cblaster gne session.json --output gne.csv --delimiter \",\"\n\n"
               "Save plot as a static HTML file:\n"
               "  $ cblaster gne session.json -p gne.html\n\n"
               "Cameron Gilchrist, 2020",
    )
    gne.add_argument(
        "session",
        type=lambda x: full_path(x, os.R_OK),
        help="cblaster session file"
    )
    add_gne_params_group(gne)
    add_gne_output_group(gne)


def add_extract_subparser(subparsers):
    parser = subparsers.add_parser(
        "extract",
        help="Extract hit sequences from session files",
        description="Extract information from session files",
        epilog="Example usage\n-------------\n"
               "Extract names of sequences matching a specific query:\n"
               "  $ cblaster extract session.json -q \"Query1\"\n\n"
               "Extract, download from NCBI and write to file in FASTA format:\n"
               "  $ cblaster extract session.json -q \"Query1\" -d -o output.fasta\n\n"
               "Extract only from specific organisms (regular expressions):\n"
               "  $ cblaster extract session.json -or \"Aspergillus.*\" \"Penicillium.*\"\n\n"
               "Generate delimited table (CSV) of all hits in clusters:\n"
               "  $ cblaster extract session.json -de \",\"\n\n"
               "Cameron Gilchrist, 2020",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("session", help="cblaster session file")

    fil = parser.add_argument_group("Filters")
    fil.add_argument(
        "-q",
        "--queries",
        help="IDs of query sequences",
        nargs="+"
    )
    fil.add_argument(
        "-or",
        "--organisms",
        help="Organism names, accepts regular expressions",
        nargs="+"
    )
    fil.add_argument(
        "-sc",
        "--scaffolds",
        help="Scaffold names/ranges, in the form scaffold_name:start-stop",
        nargs="+"
    )

    out = parser.add_argument_group("Output")
    out.add_argument(
        "-o",
        "--output",
        help="Output file name"
    )
    out.add_argument(
        "-es",
        "--extract_sequences",
        help="Extract protein sequences for all extracted proteins. The resulting summary will"
             "have a fasta format.",
        action="store_true",
    )
    out.add_argument(
        "-no",
        "--name_only",
        help="Do not save sequence descriptions (i.e. no genomic coordinates)",
        action="store_true",
    )
    out.add_argument(
        "-de",
        "--delimiter",
        help="Sequence description delimiter"
    )


def add_extract_clusters_subparser(subparsers):
    parser = subparsers.add_parser(
        "extract_clusters",
        help="Extract clusters from a session file in genbank format",
        description="Extract clusters from a session file",
        epilog="Example usage\n-------------\n"
               "Extract all clusters (this can take a while for remote sessions):\n"
               "  $ cblaster extract_clusters session.json -o example_directory\n\n"
               "Extract clusters 1 through 10 + cluster 25 (numbers found in cblaster"
               " search output):\n"
               "  $ cblaster extract_clusters session.json -c 1-10 25 -o example_directory\n\n"
               "Extract only from specific organism/s (using regular expressions):\n"
               "  $ cblaster extract_clusters session.json -or \"Aspergillus.*\" \"Penicillium.*\" "
               "-o example_directory\n\n"
               "Extract clusters from a specific range on scaffold_123 + all clusters on scaffold_234:\n"
               "  $ cblaster extract_clusters session.json -sc scaffold_123:1-80000 scaffold_234 -o "
               "example_directory\n\n"
               "Cameron Gilchrist, 2021",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "session",
        type=lambda x: full_path(x, os.R_OK),
        help="cblaster session file"
    )

    filter_parser = parser.add_argument_group("Filters")
    filter_parser.add_argument(
        "-c",
        "--clusters",
        help="Cluster numbers/ ranges provided by the summary file of the 'search' command.",
        nargs="+"
    )
    filter_parser.add_argument(
        "-st",
        "--score_threshold",
        help="Minimum score of a cluster in order to be included",
        type=float
    )
    filter_parser.add_argument(
        "-or",
        "--organisms",
        help="Organism names (can be regex patterns)",
        nargs="+"
    )
    filter_parser.add_argument(
        "-sc",
        "--scaffolds",
        help="Scaffold names/ranges e.g name:start-stop. Only clusters fully within the range are selected.",
        nargs="+"
    )

    output = parser.add_argument_group("Output options")
    output.add_argument(
        "-o",
        "--output",
        type=lambda x: full_path(x, os.W_OK, dir=True),
        help="Output directory for the clusters",
        required=True
    )
    output.add_argument(
        "-pf",
        "--prefix",
        help="Start of the name for each cluster file, the base name is 'cluster_clutser_number' e.g. cluster1",
        default=""
    )
    output.add_argument(
        "-f",
        "--format",
        choices=["genbank", "bigscape"],
        help="The format of the resulting files. The options are genbank and bigscape",
        default="genbank"
    )
    output.add_argument(
        "-mc",
        "--maximum_clusters",
        type=int,
        default=50,
        help="The maximum amount of clusters that will be extracted. Ordered on score (def. 50)"
    )


def add_plot_clusters_subparser(subparsers):
    desc = "Plot clusters using clinker"
    parser = subparsers.add_parser(
        "plot_clusters",
        help=desc,
        description=desc,
        epilog="Example usage\n-------------\n"
               "Plot all clusters (up to --maximum_clusters):\n"
               " $ cblaster plot_clusters session.json\n\n"
               "Plot clusters 1 through 10 + cluster 25 (numbers found in cblaster"
               " search output):\n"
               " $ cblaster plot_clusters session.json -c 1-10 25 -o plot.html\n\n"
               "Plot only specific organism/s (using regular expressions):\n"
               " $ cblaster plot_clusters session.json -or \"Aspergillus.*\" \"Penicillium.*\" -o plot.html\n\n"
               "Plot clusters from a specific range on scaffold_123 + all clusters on scaffold_234:\n"
               "  $ cblaster plot_clusters session.json -sc scaffold_123:1-80000 scaffold_234 -o plot.html\n\n"
               "Cameron Gilchrist, 2020",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "session",
        type=lambda x: full_path(x, os.R_OK),
        help="cblaster session file"
    )

    filter_parser = parser.add_argument_group("Filters")
    filter_parser.add_argument(
        "-c",
        "--clusters",
        help="Cluster numbers/ ranges provided by the summary file of the 'search' command.",
        nargs="+"
    )
    filter_parser.add_argument(
        "-st",
        "--score_threshold",
        help="Minimum score of a cluster to be included",
        type=float
    )
    filter_parser.add_argument(
        "-or",
        "--organisms",
        help="Organism names",
        nargs="+"
    )
    filter_parser.add_argument(
        "-sc",
        "--scaffolds",
        help="Scaffold names/ranges",
        nargs="+"
    )

    output = parser.add_argument_group("Output options")
    output.add_argument(
        "-o",
        "--output",
        type=lambda x: full_path(x, os.W_OK),
        help="Location were to store the plot file."
    )
    output.add_argument(
        "-mc",
        "--maximum_clusters",
        type=int,
        default=50,
        help="The maximum amount of clusters that will be plotted. Ordered on score (def. 50)"
    )


def add_config_subparser(subparsers):
    desc = "Configure cblaster"
    parser = subparsers.add_parser(
        "config",
        help=desc,
        description="Configure cblaster (e.g. for setting NCBI e-mail addresses or API keys)",
        epilog=(
            "Example usage\n-------------\n"
            "Set an e-mail address:\n"
            " $ cblaster config --email \"foo@bar.com\"\n\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--email",
        help="Your e-mail address, required by NCBI to prevent abuse",
        type=str,
    )
    parser.add_argument(
        "--api_key",
        help="NCBI API key",
        type=str,
    )
    parser.add_argument(
        "--max_tries",
        help="How many times failed requests are retried (def. 3)",
        type=int,
    )


def get_parser():
    parser = argparse.ArgumentParser(
        "cblaster",
        description="cblaster finds co-located sequence homologues.\n"
                    "Documentation available at https://cblaster.readthedocs.io\n"
                    "Type -h/--help after either subcommand for full description of"
                    " available arguments.",
        epilog="Cameron Gilchrist, 2020",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version", version="%(prog)s " + __version__)
    parser.add_argument("-d", "--debug", help="Print debugging information", action="store_true")
    parser.add_argument(
        "--testing",
        action="store_true",
        help="Turn this argument on when running tests to suppress certain actions like opening plots")
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
    add_extract_subparser(subparsers)
    add_extract_clusters_subparser(subparsers)
    add_plot_clusters_subparser(subparsers)
    add_config_subparser(subparsers)
    return parser


def parse_args(args):
    parser = get_parser()
    arguments = parser.parse_args(args)

    if not arguments.subcommand:
        parser.print_help()
        raise SystemExit

    if arguments.subcommand in (
        "gui",
        "makedb",
        "gne",
        "extract",
        "extract_clusters",
        "plot_clusters",
        "config",
    ):
        return arguments

    if getattr(arguments, "query_file") and getattr(arguments, "query_ids"):
        parser.error("'--query_file' and '--query_ids' are mutually exclusive")

    if arguments.mode == "remote":
        if arguments.database[0] not in NCBI_DATABASES:
            parser.error(f"Valid databases are: {', '.join(NCBI_DATABASES)}")
    else:
        for arg in ["entrez_query", "rid"]:
            if getattr(arguments, arg):
                parser.error(f"--{arg} can only be used when --mode is 'remote'")

    # Convert key to its corresponding builtin function
    arguments.binary_key = getattr(builtins, arguments.binary_key)

    if arguments.recompute and not arguments.session_file:
        parser.error("--recompute requires --session_file")

    return arguments
