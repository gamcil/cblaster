#!/usr/bin/env python3

"""
This module stores small helper functions.
"""

from pathlib import Path

import shutil
import requests


def get_program_path(aliases):
    """Get programs path given a list of program names."""
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


def prepare_query_ids(query_ids):
    """Handle query ID input.
    Accepts either an iterable of str IDs or path to newline delimited text file
    containing IDs.
    """
    if isinstance(query_ids, (list, tuple, set)):
        return query_ids
    if isinstance(query_ids, str):
        if not Path(query_ids).exists():
            raise FileNotFoundError("No file exists at given path")
        return [
            line.strip() for line in Path(query_ids).read_text().split("\n") if line
        ]
    raise ValueError("Expected an iterable of str, or path to query ID file")


def parse_fasta(handle):
    """Parse sequences in a FASTA file.

    Returns
    -------
    dict:
        Sequences in FASTA file keyed on their headers (i.e. > line)
    """
    sequences = {}
    for line in handle:
        try:
            line = line.decode().strip()
        except AttributeError:
            line = line.strip()
        if line.startswith(">"):
            header = line[1:]
            sequences[header] = ""
        else:
            sequences[header] += line
    return sequences


def efetch_sequences_request(headers):
    response = requests.post(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?",
        params={"db": "protein", "rettype": "fasta"},
        files={"id": ",".join(headers)},
    )
    if response.status_code != 200:
        raise requests.HTTPError(
            f"Error fetching sequences from NCBI [code {response.status_code}]."
            " Bad query IDs?"
        )
    return response


def efetch_sequences(headers):
    """Retrieve protein sequences from NCBI for supplied accessions.

    This function uses EFetch from the NCBI E-utilities to retrieve the sequences for
    all synthases specified in `headers`. It then calls `fasta.parse` to parse the
    returned response; note that extra processing has to occur because the returned
    FASTA will contain a full sequence description in the header line after the
    accession.

    Parameters
    ----------
    headers : list, tuple
        A list of valid NCBI sequence identifiers (accession, GI, etc). Should
        correspond to an entry in the Protein database.
    """
    response = efetch_sequences_request(headers)
    sequences = {}
    for key, value in parse_fasta(response.text.split("\n")).items():
        for header in headers:
            if header not in sequences and header in key:
                sequences[header] = value
                break
    return sequences


def get_sequences(query_file=None, query_ids=None):
    """Convenience function to get dictionary of query sequences from file or IDs."""
    if query_file and not query_ids:
        with open(query_file) as query:
            sequences = parse_fasta(query)
    elif query_ids:
        sequences = efetch_sequences(query_ids)
    else:
        raise ValueError("Expected 'query_file' or 'query_ids'")
    return sequences
