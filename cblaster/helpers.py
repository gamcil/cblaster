#!/usr/bin/env python3


import shutil
import requests
import logging

from pathlib import Path
from collections import OrderedDict


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


def parse_fasta(handle):
    """Parse sequences in a FASTA file.

    Sequence headers are trimmed after the first whitespace.

    Returns:
        Sequences in FASTA file keyed on their headers (i.e. > line)
    """
    sequences = OrderedDict()
    skip = False

    for line in handle:
        try:
            line = line.decode().strip()
        except AttributeError:
            line = line.strip()
        if line.startswith(">"):
            header = line[1:].split(" ", 1)[0]
            skip = header in sequences
            if skip:
                LOG.warning("Skipping duplicate sequence: %s", header)
            else:
                sequences[header] = ""
        else:
            if not skip:
                sequences[header] += line
    return sequences


def parse_fasta_file(path):
    with open(path) as fp:
        sequences = parse_fasta(fp)
    return sequences


def efetch_sequences_request(headers):
    """Launch E-Fetch request for a list of sequence accessions.

    Parameters:
        headers (list): NCBI sequence accessions.
    Raises:
        requests.HTTPError: Received bad status code from NCBI.
    Returns:
        requests.models.Response: Response returned by requests library.
    """
    response = requests.post(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?",
        params={"db": "protein", "rettype": "fasta"},
        files={"id": ",".join(headers)},
    )

    LOG.debug("Efetch IDs: %s", headers)
    LOG.debug("Efetch URL: %s", response.url)

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

    Parameters:
        headers (list): Valid NCBI sequence identifiers (accession, GI, etc.).
    """
    LOG.info("Querying NCBI for sequences of: %s", headers)
    response = efetch_sequences_request(headers)
    results = response.text.split("\n")
    return parse_fasta(results)


def sequences_to_fasta(sequences):
    """Formats sequence dictionary as FASTA."""
    return "\n".join(
        f">{header}\n{sequence}"
        for header, sequence in sequences.items()
    )


def get_sequences(query_file=None, query_ids=None):
    """Convenience function to get dictionary of query sequences from file or IDs.

    Parameters:
        query_file (str): Path to FASTA file containing query protein sequences.
        query_ids (list): NCBI sequence accessions.
    Raises:
        ValueError: Did not receive values for query_file or query_ids.
    Returns:
        sequences (dict): Dictionary of query sequences keyed on accession.
    """
    if query_file and not query_ids:
        with open(query_file) as query:
            sequences = parse_fasta(query)
    elif query_ids:
        sequences = efetch_sequences(query_ids)
    else:
        raise ValueError("Expected 'query_file' or 'query_ids'")
    return sequences


def get_project_root():
    return Path(__file__).resolve().parent
