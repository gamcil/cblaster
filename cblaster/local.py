#!/usr/bin/env python3

"""
This module handles BLAST search using local tools (DIAMOND and BLASTp).
"""

import logging
import subprocess

from tempfile import NamedTemporaryFile as NTF

from cblaster import helpers
from cblaster.classes import Hit


LOG = logging.getLogger(__name__)


def parse(results, min_identity=30, min_coverage=50, max_evalue=0.01):
    """Parse a string containing results of a BLAST/DIAMOND search.

    Parameters
    ----------
    results: str
        Result of `diamond()` or `blastp()`, the stdout from diamond/blastp
        subprocess call.
    min_identity: float
    min_coverage: float
    max_evalue: float

    Return
    ------
    hits: list
        Hit objects representing each hit that surpasses the scoring thresholds.
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

    Parameters
    ----------
    fasta: str
        Path to FASTA query file.
    database: str
        Path to a DIAMOND database to be searched against. This should consist of
        sequences derived from NCBI, such that `sseqid` field of results will
        correspond to valid NCBI identifiers.
    max_evalue: float
        Maximum e-value.
    min_identity: float
        Minimum percent identity.
    min_coverage: float
        Minimum percent query coverage.
    cpus: int
        Number of CPU threads to use.

    Returns
    -------
    str
        Search results, taken directly from stdout of the subprocess call.
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


def search(database, query_file=None, query_ids=None, **kwargs):
    """Launch a new BLAST search using either DIAMOND or command-line BLASTp (remote)."""
    if query_file and not query_ids:
        return _search_file(query_file, database, **kwargs)
    if query_ids:
        return _search_ids(query_ids, database, **kwargs)
    raise ValueError("Expected either 'query_ids' or 'query_file'")
