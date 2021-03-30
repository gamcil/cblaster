#!/usr/bin/env python3

import shutil
import requests
import logging
import time

from collections import OrderedDict
from pathlib import Path

from Bio import SeqIO, Entrez

from cblaster import config, genome_parsers as gp


LOG = logging.getLogger(__name__)


def get_program_path(aliases):
    """Get programs path given a list of program names.

    Parameters:
        aliases (list): Program aliases, e.g. ["diamond", "diamond-aligner"]
    Raises:
        ValueError: Could not find any of the given aliases on system $PATH.
    Returns:
        Path to program executable.
    """
    for alias in aliases:
        which = shutil.which(alias)
        if which:
            return which
    raise ValueError(f"Failed to find {aliases} on system $PATH!")


def form_command(parameters):
    """Flatten a dictionary to create a command list for use in subprocess.run()"""
    command = [] if "args" not in parameters else parameters.pop("args")
    for key, value in parameters.items():
        if isinstance(value, list):
            command.extend([key, *value])
        else:
            command.extend([key, value])
    return command


def efetch_sequences(headers):
    """Retrieve protein sequences from NCBI for supplied accessions.

    This function uses EFetch from the NCBI E-utilities to retrieve the sequences for
    all synthases specified in `headers`. The calls to EFetch can not exceed 500 accessions
    this means that the calls have to be limited. It then calls `fasta.parse` to parse the
    returned response; note that extra processing has to occur because the returned
    FASTA will contain a full sequence description in the header line after the
    accession.

    Args:
        headers (list): Valid NCBI sequence identifiers (accession, GI, etc.).
    Returns:
        a dictionary of sequences keyed on header id
    """
    LOG.info("Querying NCBI for %d sequences.", len(headers))
    try:
        handle = Entrez.efetch(
            db="protein",
            id=headers,
            rettype="fasta",
            retmode="text",
        )
    except IOError:
        LOG.exception("Failed to fetch sequences")
        raise
    return {
        record.name: str(record.seq)
        for record in SeqIO.parse(handle, 'fasta')
    }


def sequences_to_fasta(sequences):
    """Formats sequence dictionary as FASTA."""
    return "\n".join(
        f">{header}\n{sequence}"
        for header, sequence in sequences.items()
    )


def get_project_root():
    return Path(__file__).resolve().parent


def get_sequences(query_file=None, query_ids=None, query_profiles=None):
    """Convenience function to get dictionary of query sequences from file or IDs.

    Parameters:
        query_file (str): Path to FASTA genbank or EMBL file containing query
        protein sequences.
        query_ids (list): NCBI sequence accessions.
        query_profiles (list): Pfam profile accessions.
    Raises:
        ValueError: Did not receive values for query_file or query_ids.
    Returns:
        sequences (dict): Dictionary of query sequences keyed on accession.
    """
    if query_file and not query_ids:
        organism = gp.parse_file(query_file)
        if Path(query_file).suffix.lower() in gp.FASTA_SUFFIXES:
            sequences = OrderedDict((r.id, str(r.seq)) for r in organism["records"])
        else:
            genes = gp.organisms_to_tuples([organism])
            sequences = OrderedDict(
                (gene[1], gene[5])
                for gene in genes if gene[0] != "scaffold"
            )
    elif query_ids:
        sequences = efetch_sequences(query_ids)
    elif query_profiles:
        sequences = None
    else:
        raise ValueError("Expected 'query_file' or 'query_ids', or 'query_profiles'")
    return sequences
