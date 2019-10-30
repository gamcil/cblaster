#!/usr/bin/env python3

"""
This module handles the querying of NCBI for the genomic context of sequences.
"""

import logging

from collections import defaultdict

import requests
import operator

from clusterblaster.classes import Organism, Scaffold


LOG = logging.getLogger(__name__)


def efetch_IPGs(ids, output_handle=None):
    """Query NCBI for Identical Protein Groups (IPG) of query sequences.

    db = protein
    id = ',' joined list
    format = ipg
    """

    response = requests.post(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?",
        params={"db": "protein", "rettype": "ipg", "retmode": "text"},
        data={"id": ",".join(ids)},
    )

    LOG.debug("IPG search URL: %s", response.url)
    LOG.debug("IPG search IDs: %s", ids)

    if response.status_code != 200:
        raise requests.HTTPError(
            "Error fetching sequences from NCBI [code {response.status_code}]."
        )

    if output_handle:
        LOG.info("Writing IPG output to %s", output_handle.name)
        output_handle.write(response.text)

    return response


def parse_IPG_table(results_handle, hits):
    """Parse the results of an IPG query.

    Parameters
    ----------
    results_handle : open file handle
        Results from identical protein groups search.
    hits : list
        Hits that were used to query NCBI.
    conserve: int
    gap: int

    Returns
    -------
    organisms : list
        Organism objects containing hit information.
    """

    hits_dict = {hit.subject: hit for hit in hits}

    organisms = defaultdict(dict)
    _ipg = None

    for line in results_handle:
        if not line or line.startswith("Id\tSource") or line.isspace():
            continue

        ipg, _, accession, start, end, strand, protein, _, organism, strain, assembly = line.split(
            "\t"
        )

        if not accession or not assembly:
            # Avoid vectors, single Gene nucleotide entries, etc
            continue

        if ipg == _ipg:
            continue
        _ipg = ipg

        # Some org. names contain strain, try to remove
        if strain in organism:
            organism = organism.replace(strain, "").strip()

        # Only make new Organism instance if not already one of this name/strain
        if strain not in organisms[organism]:
            LOG.debug("New organism: %s %s", organism, strain)
            organisms[organism][strain] = Organism(name=organism, strain=strain)

        # Add new scaffold
        if accession not in organisms[organism][strain].scaffolds:
            LOG.debug("New scaffold: %s", accession)
            organisms[organism][strain].scaffolds[accession] = Scaffold(accession)

        try:
            hit = hits_dict.pop(protein)
        except KeyError:
            continue

        # Update hit with genomic position
        hit.end = int(end)
        hit.start = int(start)
        hit.strand = strand

        # Add to corresponding scaffold
        organisms[organism][strain].scaffolds[accession].hits.append(hit)

    return [organism for strains in organisms.values() for organism in strains.values()]


def find_clusters_in_hits(hits, conserve=3, gap=20000):
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
    if conserve < 0 or gap < 0:
        raise ValueError("Expected positive integer")

    total_hits = len(hits)

    if total_hits < conserve:
        return []

    hits.sort(key=operator.attrgetter("start"))

    i, groups = 0, []
    while i < total_hits:
        for j in range(i + 1, total_hits):
            grouped = hits[j].start - hits[j - 1].end <= gap

            if j == total_hits - 1:
                # last element
                # test if a) part of previous group, or b) alone
                # then,   a) add previous group W/ last element
                #         b) add previous + last as separate
                if grouped:
                    groups.append(hits[i:])
                else:
                    groups.extend([hits[i:j], hits[j:]])

                i = total_hits
                break

            if not grouped:
                groups.append(hits[i:j])
                i = j
                break

    # Only save groups that have enough conserved QUERY hits
    # i.e. filter out blocks of hits only related to one query sequence
    return [group for group in groups if len(set(h.query for h in group)) >= conserve]


def find_clusters_in_organism(organism, conserve=3, gap=20000):
    """Run find_clusters() on all Hits on Scaffolds in an Organism instance."""
    for scaffold in organism.scaffolds.values():
        scaffold.clusters = find_clusters_in_hits(scaffold.hits, conserve, gap)

        LOG.debug(
            "Organism: %s, Scaffold: %s, Clusters: %i",
            organism.full_name,
            scaffold.accession,
            len(scaffold.clusters),
        )


def search(hits, conserve, gap):
    """Get genomic context for a collection of BLAST hits."""
    groups = efetch_IPGs([hit.subject for hit in hits])

    organisms = parse_IPG_table(groups.text.split("\n"), hits)

    LOG.info("Finding hit clusters in %i organisms", len(organisms))
    for organism in organisms:
        find_clusters_in_organism(organism, conserve, gap)

    return organisms
