#!/usr/bin/env python3


import logging
import subprocess
import os

from tempfile import NamedTemporaryFile as NTF
from pathlib import Path

from cblaster import helpers
from cblaster.classes import Hit


LOG = logging.getLogger(__name__)


def parse(results, min_identity=30, min_coverage=50, max_evalue=0.01):
    """Parse a string containing results of a BLAST/DIAMOND search.

    Arguments:
        results (list): Results returned by diamond() or blastp()
        min_identity (float): Minimum identity (%) cutoff
        min_coverage (float): Minimum coverage (%) cutoff
        max_evalue (float): Maximum e-value threshold
    Returns:
        list: Hit objects representing hits that surpass scoring thresholds
    """
    hits = []
    for row in results[:-1]:
        hit = Hit(*row.split("\t"))
        if (
            hit.identity > min_identity
            and hit.coverage > min_coverage
            and hit.evalue < max_evalue
        ):
            hits.append(hit)
    if len(hits) == 0:
        raise SystemExit("No results found")
    return hits


def diamond(fasta, database, max_evalue=0.01, min_identity=30, min_coverage=50, cpus=1):
    """Launch a local DIAMOND search against a database.

    Arguments:
        fasta (str): Path to FASTA format query file
        database (str): Path to DIAMOND database generated with cblaster makedb
        max_evalue (float): Maximum e-value threshold
        min_identity (float): Minimum identity (%) cutoff
        min_coverage (float): Minimum coverage (%) cutoff
        cpus (int): Number of CPU threads for DIAMOND to use
    Returns:
        list: Rows from DIAMOND search result table (split by newline)
    """
    diamond = helpers.get_program_path(["diamond", "diamond-aligner"])
    LOG.debug("diamond path: %s", diamond)

    parameters = {
        "args": [diamond, "blastp"],
        "--query": fasta,
        "--db": database,
        "--id": str(min_identity),
        "--evalue": str(max_evalue),
        "--outfmt": [
            "6",
            "qseqid",
            "sseqid",
            "pident",
            "qcovhsp",
            "evalue",
            "bitscore",
        ],
        "--threads": str(cpus),
        "--query-cover": str(min_coverage),
        "--max-hsps": "1",
    }

    command = helpers.form_command(parameters)
    LOG.debug("Parameters: %s", command)

    results = subprocess.run(
        command,
        stderr=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        check=True
    )

    return results.stdout.decode().split("\n")


def search(
    database,
    sequences=None,
    query_file=None,
    query_ids=None,
    blast_file=None,
    **kwargs,
):
    """Launch a new BLAST search using either DIAMOND or command-line BLASTp (remote).

    Arguments:
        database (str): Path to DIAMOND database
        query_file (str): Path to FASTA file containing query sequences
        query_ids (list): NCBI sequence accessions
    Raises:
        ValueError: No value given for query_file or query_ids
    Returns:
        list: Parsed rows with hits from DIAMOND results table
    """
    if query_file:
        table = diamond(query_file, database, **kwargs)
    else:
        if not sequences:
            sequences = helpers.get_sequences(query_ids=query_ids)

        # delete=False since you cannot open tempfiles twice in Windows
        fasta = NTF("w", delete=False)
        text = helpers.sequences_to_fasta(sequences)
        try:
            with fasta:
                fasta.write(text)
            table = diamond(fasta.name, database, **kwargs)
        finally:
            os.unlink(fasta.name)

    results = parse(table)

    if blast_file:
        LOG.info("Writing DIAMOND hit table to %s", blast_file.name)
        blast = "\n".join(results)
        blast_file.write(blast)

    return results
