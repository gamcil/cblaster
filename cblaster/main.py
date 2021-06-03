#!/usr/bin/env python3


import logging
import sys

from pathlib import Path

from Bio import Entrez

from cblaster import (
    config,
    context,
    database,
    extract,
    extract_clusters,
    helpers,
    hmm_search,
    local,
    parsers,
    plot_clusters,
    remote,
)
from cblaster.classes import Session
from cblaster.plot import plot_session, plot_gne
from cblaster.formatters import summarise_gne
from cblaster.intermediate_genes import find_intermediate_genes


logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s - %(message)s",
    datefmt="%H:%M:%S",
)
# make sure to not configure a name otherwise a different logger instance is returned where the debug level is set
# resulting in no debug information being printed
LOG = logging.getLogger()


def set_entrez():
    """Set the Entrez parameters from config"""
    cfg = config.get_config_parser()
    Entrez.email = cfg["cblaster"].get("email", None)
    Entrez.api_key = cfg["cblaster"].get("api_key", None)


def gne(
    session,
    output=None,
    max_gap=100000,
    samples=100,
    scale="linear",
    plot=None,
    hide_headers=False,
    delimiter=",",
    decimals=4,
    testing=False,
):
    """Estimate gene neighbourhood."""
    LOG.info("Starting cblaster gene neighbourhood estimation")
    LOG.info("Loading session from: %s", session)
    with open(session) as fp:
        session = Session.from_json(fp)

    LOG.info("Computing gene neighbourhood statistics")
    results = context.estimate_neighbourhood(
        session,
        max_gap=max_gap,
        samples=samples,
        scale=scale
    )
    if output:
        LOG.info("Writing GNE table to %s", output)
        summary = summarise_gne(
            results,
            hide_headers=hide_headers,
            delimiter=delimiter,
            decimals=decimals,
        )
        with open(output, "w") as f:
            f.write(summary)

    plot_gne(results, output=plot, testing=testing)
    LOG.info("Done.")


def cblaster(
    query_file=None,
    query_ids=None,
    query_profiles=None,
    mode=None,
    databases=None,
    database_pfam=None,
    gap=20000,
    unique=3,
    min_hits=3,
    min_identity=30,
    min_coverage=50,
    max_evalue=0.01,
    percentage=None,
    entrez_query=None,
    output=None,
    output_hide_headers=False,
    output_delimiter=None,
    output_decimals=4,
    output_sort_clusters=False,
    binary=None,
    binary_hide_headers=True,
    binary_delimiter=None,
    binary_key=len,
    binary_attr="identity",
    binary_decimals=4,
    rid=None,
    require=None,
    session_file=None,
    indent=None,
    plot=False,
    max_plot_clusters=50,
    recompute=False,
    blast_file=None,
    ipg_file=None,
    hitlist_size=None,
    cpus=None,
    intermediate_genes=False,
    intermediate_gene_distance=5000,
    intermediate_max_clusters=100,
    testing=False,
):
    """Run cblaster.

    This function is the central workflow for the entire cblaster package.

    Arguments:
        query_file (str): Path to FASTA format query file
        query_ids (list): NCBI protein sequence identifiers
        query_profiles(list): Pfam profile identifiers
        mode (str): Search mode ('local' or 'remote')
        databases (str): Search database (NCBI if remote, DIAMOND if local)
        database_pfam (str): Path to pfam db or where to download it
        gap (int): Maximum gap (kilobase) between cluster hits
        unique (int): Minimum number of query sequences with hits in clusters
        min_hits (int): Minimum number of hits in clusters
        min_identity (float): Minumum identity (%) cutoff
        min_coverage (float): Minumum coverage (%) cutoff
        max_evalue (float): Maximum e-value threshold
        percentage (int): % of query genes needed to be present in cluster
        entrez_query (str): NCBI Entrez query to filter search database
        output (str): Path to cblaster summary output file
        output_hide_headers (bool): Hide headers in summary table
        output_delimiter (str): Delimiter used in summary table
        output_decimals (int): Total decimal places in hit scores in summary table
        output_sort_clusters (bool): If the clusters in the final summary table need to sorted
        binary (str): Path to cblaster binary output file
        binary_hide_headers (bool): Hide headers in binary table
        binary_delimiter (str): Delimiter used in binary table
        binary_key (str): Key function used in binary table (len, max or sum)
        binary_attr (str): Hit attribute used for calculating cell values in binary table
        binary_decimals (int): Total decimal places in cell values in binary table
        rid (str): NCBI BLAST search request identifier (RID)
        require (list): Query sequences that must be in hit clusters
        session_file (str): Path to cblaster session JSON file
        indent (int): Total spaces to indent JSON files
        plot (str): Path to cblaster plot HTML file
        max_plot_clusters (int): maximum clusters that are plotted when -osc (sort on score ) argument is used
        recompute (str): Path to recomputed session JSON file
        blast_file (str): path to file to save blast output
        ipg_file (str): path to file to save ipg output
        cpus (int): number of cpu's to use when blasting.
        intermediate_genes (bool): Signifies if intermediate genes have to be shown
        hitlist_size (int): Number of database sequences to keep
        intermediate_gene_distance (int): the maximum allowed distance between the
         edge of a cluster and an intermediate gene.
        intermediate_max_clusters (int): the maximum amount of clusters for which intermediate
         genes will be fetched, since this can become expensive for remote searches
        testing (bool): flag to make sure certain code does not run when testing

    Returns:
        Session: cblaster search Session object
    """
    if session_file and all(Path(sf).exists() for sf in session_file):
        LOG.info("Loading session(s) %s", session_file)
        session = Session.from_files(session_file)

        if recompute:
            LOG.info("Filtering session with new thresholds")
            context.filter_session(
                session,
                min_identity,
                min_coverage,
                max_evalue,
                gap,
                unique,
                min_hits,
                require,
                percentage,
            )

            if intermediate_genes:
                find_intermediate_genes(
                    session, intermediate_gene_distance, intermediate_max_clusters
                )

            if recompute is not True:
                LOG.info("Writing recomputed session to %s", recompute)
                session.params["min_identity"] = min_identity
                session.params["min_coverage"] = min_coverage
                session.params["max_evalue"] = max_evalue
                session.params["require"] = require
                with open(recompute, "w") as fp:
                    session.to_json(fp, indent=indent)
    else:
        # Create a cblaster Cluster object from query input
        query = helpers.parse_query_sequences(
            query_file=query_file,
            query_ids=query_ids,
            query_profiles=query_profiles,
        )
        
        # Create a cblaster Session
        session = Session(
            query=query,
            queries=query.names,
            params={
                "mode": mode,
                "database": databases,
                "min_identity": min_identity,
                "min_coverage": min_coverage,
                "max_evalue": max_evalue,
                "require": require,
            },
        )

        if query_file:
            # get_sequences() returns OrderedDict, so save keys to
            # preserve query order
            session.params["query_file"] = query_file

        sqlite_db = None
        session.params["rid"] = rid

        if "combi" in mode and not len(databases) == 2:
            raise RuntimeError("Expected two databases for 'combi_' modes")

        if mode in ("hmm", "combi_local", "combi_remote"):
            sqlite_db = helpers.find_sqlite_db(databases[0])
            results = hmm_search.perform_hmmer(
                fasta=databases[0],
                query_profiles=query_profiles,
                pfam=database_pfam,
                session=session
            )

            # Delete first (FASTA) database when doing combined searches
            # Expect .dmnd/NCBI database name for local/remote, respectively
            if "combi" in mode:
                del databases[0]

            LOG.info("Found %i hits meeting score thresholds for hmm search", len(results))
            LOG.info("Fetching genomic context of hits")
            organisms = context.search(
                results,
                sqlite_db=sqlite_db,
                unique=unique,
                min_hits=min_hits,
                gap=gap,
                require=require,
                ipg_file=ipg_file,
                query_sequence_order=session.queries,
                percentage=percentage,
            )
            session.organisms.extend(organisms)

        # When running combined modes, run local/remote search right after HMM search
        if mode == "combi_local":
            mode = "local"
        elif mode == "combi_remote":
            mode = "remote"

        if mode == "local":
            LOG.info("Starting cblaster in local mode")
            sqlite_db = helpers.find_sqlite_db(databases[0])
            results = local.search(
                databases[0],
                sequences=session.query.sequences,
                min_identity=min_identity,
                min_coverage=min_coverage,
                max_evalue=max_evalue,
                blast_file=blast_file,
                cpus=cpus,
            )
            LOG.info(
                "Found %i hits meeting score thresholds for local search", len(results)
            )
            LOG.info("Fetching genomic context of hits")
            organisms = context.search(
                results,
                sqlite_db=sqlite_db,
                unique=unique,
                min_hits=min_hits,
                gap=gap,
                require=require,
                ipg_file=ipg_file,
                query_sequence_order=session.queries,
                percentage=percentage,
            )
            session.organisms.extend(organisms)

        elif mode == "remote":
            LOG.info("Starting cblaster in remote mode")

            # Set up mandatory Entrez params
            set_entrez()

            if not Entrez.email and not Entrez.api_key:
                raise IOError("No e-mail or NCBI API key found, please run cblaster config")

            if entrez_query:
                session.params["entrez_query"] = entrez_query
            rid, results = remote.search(
                sequences=session.query.sequences,
                rid=rid,
                database=databases[0],
                min_identity=min_identity,
                min_coverage=min_coverage,
                max_evalue=max_evalue,
                entrez_query=entrez_query,
                blast_file=blast_file,
                hitlist_size=hitlist_size,
            )
            session.params["rid"] = rid
            LOG.info(
                "Found %i hits meeting score thresholds for remote search", len(results)
            )
            LOG.info("Fetching genomic context of hits")
            organisms = context.search(
                results,
                unique=unique,
                min_hits=min_hits,
                gap=gap,
                require=require,
                ipg_file=ipg_file,
                query_sequence_order=session.queries,
                percentage=percentage,
            )
            session.organisms.extend(organisms)

        if sqlite_db:
            session.params["sqlite_db"] = str(sqlite_db)

        if intermediate_genes:
            find_intermediate_genes(
                session, intermediate_gene_distance, intermediate_max_clusters
            )

        if session_file:
            LOG.info("Writing current search session to %s", session_file[0])
            if len(session_file) > 1:
                LOG.warning("Multiple session files specified, using first")
            with open(session_file[0], "w") as fp:
                session.to_json(fp, indent=indent)

    if binary:
        LOG.info("Writing binary summary table to %s", binary)
        session.format(
            "binary",
            open(binary, "w"),
            hide_headers=binary_hide_headers,
            delimiter=binary_delimiter,
            key=binary_key,
            attr=binary_attr,
            decimals=binary_decimals,
            sort_clusters=output_sort_clusters,
        )

    LOG.info("Writing summary to %s", "stdout" if output is None else output)
    session.format(
        "summary",
        fp=open(output, "w") if output else sys.stdout,
        hide_headers=output_hide_headers,
        delimiter=output_delimiter,
        decimals=output_decimals,
        sort_clusters=output_sort_clusters,
    )

    if plot:
        plot = None if plot is True else plot
        plot_session(
            session,
            output=plot,
            sort_clusters=output_sort_clusters,
            max_clusters=max_plot_clusters,
            testing=testing,
        )

    LOG.info("Done.")
    return session


def main():
    """cblaster entry point."""
    args = parsers.parse_args(sys.argv[1:])

    if args.debug:
        LOG.setLevel(logging.DEBUG)

    if args.subcommand == "makedb":
        database.makedb(
            args.paths,
            database=args.name,
            cpus=args.cpus,
            batch=args.batch,
            force=args.force,
        )

    elif args.subcommand == "search":
        cblaster(
            query_file=args.query_file,
            query_ids=args.query_ids,
            query_profiles=args.query_profiles,
            mode=args.mode,
            databases=args.database,
            database_pfam=args.database_pfam,
            gap=args.gap,
            unique=args.unique,
            min_hits=args.min_hits,
            require=args.require,
            percentage=args.percentage,
            min_identity=args.min_identity,
            min_coverage=args.min_coverage,
            max_evalue=args.max_evalue,
            entrez_query=args.entrez_query,
            output=args.output,
            output_hide_headers=args.output_hide_headers,
            output_delimiter=args.output_delimiter,
            output_decimals=args.output_decimals,
            output_sort_clusters=args.sort_clusters,
            binary=args.binary,
            binary_hide_headers=args.binary_hide_headers,
            binary_delimiter=args.binary_delimiter,
            binary_key=args.binary_key,
            binary_attr=args.binary_attr,
            binary_decimals=args.binary_decimals,
            rid=args.rid,
            session_file=args.session_file,
            indent=args.indent,
            recompute=args.recompute,
            plot=args.plot,
            max_plot_clusters=args.max_plot_clusters,
            blast_file=args.blast_file,
            ipg_file=args.ipg_file,
            hitlist_size=args.hitlist_size,
            cpus=args.cpus,
            intermediate_genes=args.intermediate_genes,
            intermediate_gene_distance=args.max_distance,
            intermediate_max_clusters=args.maximum_clusters,
            testing=args.testing,
        )

    elif args.subcommand == "gui":
        from cblaster.gui.main import cblaster_gui

        cblaster_gui()

    elif args.subcommand == "gne":
        gne(
            args.session,
            args.output,
            max_gap=args.max_gap,
            samples=args.samples,
            scale=args.scale,
            delimiter=args.delimiter,
            hide_headers=args.hide_headers,
            decimals=args.decimals,
            plot=args.plot,
            testing=args.testing,
        )

    elif args.subcommand == "extract":
        extract.extract(
            args.session,
            extract_seqs=args.extract_sequences,
            output=args.output,
            queries=args.queries,
            organisms=args.organisms,
            scaffolds=args.scaffolds,
            name_only=args.name_only,
            delimiter=args.delimiter,
        )

    elif args.subcommand == "extract_clusters":
        set_entrez()
        extract_clusters.extract_clusters(
            args.session,
            args.output,
            prefix=args.prefix,
            cluster_numbers=args.clusters,
            score_threshold=args.score_threshold,
            organisms=args.organisms,
            scaffolds=args.scaffolds,
            format_=args.format,
            max_clusters=args.maximum_clusters,
        )

    elif args.subcommand == "plot_clusters":
        set_entrez()
        plot_clusters.plot_clusters(
            session=args.session,
            cluster_numbers=args.clusters,
            score_threshold=args.score_threshold,
            organisms=args.organisms,
            scaffolds=args.scaffolds,
            plot_outfile=args.output,
            max_clusters=args.maximum_clusters,
            testing=args.testing,
        )

    elif args.subcommand == "config":
        if not args.email and not args.api_key:
            LOG.info(
                "No e-mail or API key specified; if this is your first time"
                " running cblaster config, please make sure you provide one."
            )
        config.write_config_file(
            email=args.email,
            api_key=args.api_key,
            max_tries=args.max_tries,
        )


if __name__ == "__main__":
    main()
