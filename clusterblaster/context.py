#!/usr/bin/env python3

"""
This module handles the querying of NCBI for the genomic context of sequences.
"""

import requests
import operator

from clusterblaster.classes import Organism, Scaffold


def find_clusters(hits, conserve=3, gap=20000):
    """Find conserved hit blocks in a collection of BLAST hits.

    Parameters
    ----------
    hits: list
        Hit objects with positional information.
    conserve: int
        Minimum number of hits that must be in one cluster to be saved.
    gap: int
        Maximum intergenic distance (bp) between any two hits. This is calculated
        from the end of one gene to the start of the next. If this distance exceeds the
        specified value, the group is considered finished and checked against the
        conserve value.

    Returns
    -------
    groups: list
        Groups of Hit objects.
    """

    total_hits = len(hits)
    if total_hits < conserve:
        return []

    hits.sort(key=operator.attrgetter("start"))

    i, groups = 0, []
    while i < total_hits:
        reset, last = False, i
        for j in range(i + 1, total_hits):
            if hits[j].end - hits[j - 1].start > gap:
                reset, last = True, j - 1  # Gap is too big
            elif j == total_hits - 1:
                reset, last = True, j  # At the last element
            if reset:
                if last + 1 - i >= conserve:
                    groups.append(hits[i : last + 1])  # Save indices
                break
        i = last + 1  # Move index ahead of last block

    return groups


def efetch_IPGs(hits, output_handle=None):
    """Query NCBI for Identical Protein Groups (IPG) of query sequences.

    db = protein
    id = ',' joined list
    format = ipg
    """

    response = requests.post(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?",
        params={"db": "protein", "rettype": "ipg", "retmode": "text"},
        data={"id": ",".join(hit.subject for hit in hits)},
    )

    if response.status_code != 200:
        raise requests.HTTPError(
            "Error fetching sequences from NCBI [code {response.status_code}]."
        )

    if output_handle:
        output_handle.write(response.text)

    return response


def parse_IPG_table(results_handle, hits, conserve=3, gap=20000):
    """Parse the results of an IPG query.

    Parameters
    ----------
    results_handle : open file handle
        Results from identical protein groups search.
    hits : list
        Hits that were used to query NCBI.

    Returns
    -------
    organisms : list
        Organism objects containing hit information.
    """

    hits_dict = {hit.subject: hit for hit in hits}

    organisms = {}
    _ipg = None

    for line in results_handle:
        if not line or line.startswith("Id\tSource") or line.isspace():
            continue

        ipg, _, accession, start, end, _, protein, _, organism, strain, assembly = line.split(
            "\t"
        )

        if not accession or not assembly:
            # Avoid vectors, single Gene nucleotide entries, etc
            continue

        if ipg == _ipg:
            continue
        _ipg = ipg

        if organism not in organisms:
            organisms[organism] = Organism(name=organism, strain=strain)

        if accession not in organisms[organism].scaffolds:
            organisms[organism].scaffolds[accession] = Scaffold(accession)

        try:
            hit = hits_dict.pop(protein)
        except KeyError:
            continue

        hit.start, hit.end = int(start), int(end)
        organisms[organism].scaffolds[accession].hits.append(hit)

    for organism in organisms.values():
        for scaffold in organism.scaffolds.values():
            scaffold.clusters = find_clusters(scaffold.hits, conserve, gap)

    return organisms


def search(hits, conserve, gap):
    """Get genomic context for a collection of BLAST hits."""
    groups = efetch_IPGs(hits)
    return parse_IPG_table(groups.text.split("\n"), hits, conserve, gap)
