#!/usr/bin/env python3


import logging
import subprocess

from tempfile import NamedTemporaryFile as NTF

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
        command, stderr=subprocess.DEVNULL, stdout=subprocess.PIPE, check=True
    )

    return results.stdout.decode().split("\n")


def _search_file(fasta, database, **kwargs):
    """Launcher function for `diamond` and `blastp` modes."""
    LOG.info("Starting DIAMOND search")
    return parse(diamond(fasta, database, **kwargs))


def _search_ids(ids, database, **kwargs):
    """Thin wrapper around _search_file to facilitate query IDs instead of FASTA.

    Since _local_BLAST takes a file path as input, can pass it the name attribute of a
    NamedTemporaryFile. So, first obtain the sequences for each ID from NCBI via
    efetch_sequences, then write those to an NTF and pass to _local_BLAST.
    """
    with NTF("w") as fasta:
        for header, sequence in helpers.efetch_sequences(ids).items():
            fasta.write(f">{header}\n{sequence}\n")
        fasta.seek(0)
        results = _search_file(fasta.name, database, **kwargs)
    return results


def search(database, query_file=None, query_ids=None, blast_file=None, **kwargs):
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
    if query_file and not query_ids:
        results = _search_file(query_file, database, **kwargs)
    elif query_ids:
        results = _search_ids(query_ids, database, **kwargs)
    else:
        raise ValueError("Expected either 'query_ids' or 'query_file'")
    if blast_file:
        LOG.info("Writing DIAMOND hit table to %s", blast_file.name)
        blast = "\n".join(results)
        blast_file.write(blast)
    return results
