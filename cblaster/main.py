#!/usr/bin/env python3


import argparse
import logging
import sys
from pathlib import Path

from cblaster import (
    __version__,
    context,
    database,
    helpers,
    local,
    remote,
    parsers,
)
from cblaster.classes import Session
from cblaster.plot import plot_session


logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s - %(message)s",
    datefmt="%H:%M:%S"
)
LOG = logging.getLogger(__name__)


def makedb(genbanks, filename, indent=None):
    """Generate JSON and diamond databases."""
    LOG.info("Starting cblaster makedb")
    db = database.Database.from_files(genbanks)

    LOG.info("Writing FASTA file with database sequences: %s", filename + ".faa")
    LOG.info("Building DIAMOND database: %s", filename + ".dmnd")
    db.makedb(filename)

    LOG.info("Building JSON database: %s", filename + ".json")
    with open(f"{filename}.json", "w") as handle:
        db.to_json(handle, indent=indent)

    LOG.info("Done.")



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
    output_hide_headers=False,
    output_delimiter=None,
    output_decimals=4,
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
    save_html=None,
    recompute=False,
):
    """Run cblaster."""

    if session_file and Path(session_file).exists():
        LOG.info("Loading %s", session_file)
        with open(session_file) as fp:
            session = Session.from_json(fp)

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
        LOG.info("Writing binary summary table to %s", binary)
        session.format(
            "binary",
            open(binary, "w"),
            hide_headers=binary_hide_headers,
            delimiter=binary_delimiter,
            key=binary_key,
            attr=binary_attr,
            decimals=binary_decimals,
        )

    LOG.info("Writing summary to %s", "stdout" if output == sys.stdout else output)
    results = session.format(
        "summary",
        fp=open(output, "w") if output else sys.stdout,
        hide_headers=output_hide_headers,
        delimiter=output_delimiter,
        decimals=output_decimals,
    )

    if plot:
        plot = None if plot is True else plot
        plot_session(session, output=plot)

    LOG.info("Done.")

    return session


def main():
    args = parsers.parse_args(sys.argv[1:])

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
            output_hide_headers=args.output_hide_headers,
            output_delimiter=args.output_delimiter,
            output_decimals=args.output_decimals,
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
        )

    elif args.subcommand == "gui":
        from cblaster.gui.main import cblaster_gui
        cblaster_gui()


if __name__ == "__main__":
    main()
