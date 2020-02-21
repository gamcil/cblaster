#!/usr/bin/env python3

"""
This module handles the querying of NCBI for the genomic context of sequences.
"""

import logging
import re

from collections import defaultdict, namedtuple
from operator import attrgetter

import requests

from cblaster import database
from cblaster.classes import Organism, Scaffold


LOG = logging.getLogger(__name__)


def efetch_IPGs(ids, output_handle=None):
    """Query NCBI for Identical Protein Groups (IPG) of query sequences.

    db = protein
    id = ',' joined list
    format = ipg
    """

    # Split into chunks since retmax=10000
    chunks = [ids[i : i + 10000] for i in range(0, len(ids), 10000)]

    table = ""
    for ix, chunk in enumerate(chunks, 1):
        response = requests.post(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?",
            params={
                "db": "protein",
                "rettype": "ipg",
                "retmode": "text",
                "retmax": 10000,
            },
            data={"id": ",".join(chunk)},
        )

        if response.status_code != 200:
            raise requests.HTTPError(
                "Error fetching sequences from NCBI [code {response.status_code}]."
            )

        table += response.text

    if output_handle:
        LOG.info("Writing IPG output to %s", output_handle.name)
        output_handle.write(table)

    return table.split("\n")


def parse_IPG_table(results_handle, hits):
    """Parse the results of an IPG query.

    For every identical protein group, this function will choose the first that:
        1) Has non-empty nuccore accession, start, end and strand fields
        2) Is from either RefSeq/GenBank (INSDC)

    As per NCBI documentation, the 'best' entry is chosen in preference of RefSeq ->
    SwissProt -> PIR, PDB -> GenBank -> Patent. However, entries not from RefSeq/GenBank
    are protein-only accessions and have no corresponding nucleotide accession.
    I don't know the exact order this table is formulated, but it tends to reflect the
    above.

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

    fields = [
        "source",
        "scaffold",
        "start",
        "end",
        "strand",
        "protein_id",
        "protein_name",
        "organism",
        "strain",
        "assembly",
    ]
    Entry = namedtuple("Entry", fields)

    seen = set()
    groups = defaultdict(list)

    # Parse table, form groups of identical proteins
    for line in results_handle:
        if not line or line.startswith("Id\tSource") or line.isspace():
            continue

        ipg, *fields = line.split("\t")
        entry = Entry(*fields)

        if not all(
            [entry.scaffold, entry.start, entry.end, entry.strand]
        ) or entry.source not in ("RefSeq", "INSDC",):
            # Avoid vectors, single Gene nucleotide entries, etc
            continue

        groups[ipg].append(entry)

    organisms = defaultdict(dict)

    for ipg in list(groups):
        group = groups.pop(ipg)  # Remove groups as we go
        hit = None

        # Find a representative Hit object per IPG
        for entry in group:
            try:
                hit = hits_dict[entry.protein_id]
            except KeyError:
                continue
            break

        if not hit:
            LOG.warning("Could not find Hit object for IPG %s", ipg)
            continue

        # Now populate the organisms dictionary with copies
        for entry in group:
            if entry.protein_id in seen:
                continue

            seen.add(entry.protein_id)

            org, st, acc = entry.organism, entry.strain, entry.scaffold

            # Some organism names contain the strain, so try to remove
            if st in org:
                org = org.replace(st, "").strip()

            # If we haven't seen this strain before, create new Organism
            if st not in organisms[org]:
                organisms[org][st] = Organism(name=org, strain=st)

            # Likewise with Scaffold
            if acc not in organisms[org][st].scaffolds:
                organisms[org][st].scaffolds[acc] = Scaffold(acc)

            # Copy the representative hit, but change the subject and add location
            h = hit.copy(
                subject=entry.protein_id,
                end=int(entry.end),
                start=int(entry.start),
                strand=entry.strand,
            )

            # Save it
            organisms[org][st].scaffolds[acc].hits.append(h)

    return [organism for strains in organisms.values() for organism in strains.values()]


def query_local_DB(hits, database):
    """Build Organisms/Scaffolds using database.DB instance.

    Retrieves corresponding Protein object from the DB by passing it the Hit subject
    attribute, which should be a 4 field '|' delimited header string containing the full
    lineage of the protein (e.g. organism|strain|scaffold|protein).

    Internally uses defaultdict to build up hierarchy of unique organisms, then returns
    list of just each organism dictionary.
    """

    organisms = defaultdict(dict)

    for hit in hits:
        protein = database.proteins[hit.subject]

        # Mostly for readability..
        org, strain, scaf = protein.organism, protein.strain, protein.scaffold

        if strain not in organisms[org]:
            organisms[org][strain] = Organism(org, strain)

        if scaf not in organisms[org][strain].scaffolds:
            organisms[org][strain].scaffolds[scaf] = Scaffold(scaf)

        # Want to report just protein ID, not lineage
        hit.subject = protein.id

        # Save genomic location on the Hit instance
        hit.start = protein.start
        hit.end = protein.end
        hit.strand = protein.strand

        organisms[protein.organism][protein.strain].scaffolds[
            protein.scaffold
        ].hits.append(hit)

    return [organism for strains in organisms.values() for organism in strains.values()]


def hits_contain_required_queries(hits, queries):
    """Check that a group of Hit objects contains a Hit for each required query."""
    bools = [False] * len(queries)
    for index, query in enumerate(queries):
        for hit in hits:
            if hit.query == query:
                bools[index] = True
                break
    return all(bools)


def hits_contain_unique_queries(hits, conserve=3):
    """Check that a group of Hit objects belong to some threshold of unique queries."""
    return len(hits) >= conserve


def find_clusters(hits, require=None, conserve=3, gap=20000):
    """Find clustered Hit objects.

    Parameters
    ----------
    hits: list
        Hit objects with positional information.
    require: iterable
        Names of query sequences that must be represented in a Hit cluster.
    conserve: int
        Minimum number of unique queries hit in a Hit cluster.
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
    if total_hits == 1:
        if conserve == 1:
            return [hits]
        return []

    def conditions_met(group):
        """Check require/conserve params are satisfied."""
        req = hits_contain_required_queries(group, require) if require else True
        con = hits_contain_unique_queries(group, conserve)
        return req and con

    sorted_hits = sorted(hits, key=attrgetter("start"))
    first = sorted_hits.pop(0)
    group, border = [first], first.end

    for hit in sorted_hits:
        if hit.start <= border + gap:
            group.append(hit)
            border = max(border, hit.end)
        else:
            if conditions_met(group):
                yield group
            group, border = [hit], hit.end

    if conditions_met(group):
        yield group


def find_clusters_in_organism(organism, conserve=3, gap=20000, require=None):
    """Run find_clusters() on all Hits on Scaffolds in an Organism instance."""
    for scaffold in organism.scaffolds.values():
        scaffold.clusters = list(
            find_clusters(scaffold.hits, conserve=conserve, gap=gap, require=require)
        )
        LOG.debug(
            "Organism: %s, Scaffold: %s, Clusters: %i",
            organism.full_name,
            scaffold.accession,
            len(scaffold.clusters),
        )


def search(hits, conserve, gap, require=None, json_db=None):
    """Get genomic context for a collection of BLAST hits."""
    if json_db:
        LOG.info("Loading JSON database: %s", json_db)
        db = database.DB.from_json(json_db)
        organisms = query_local_DB(hits, db)
    else:
        rows = efetch_IPGs([hit.subject for hit in hits])
        organisms = parse_IPG_table(rows, hits)

    LOG.info("Searching for clustered hits across %i organisms", len(organisms))
    for organism in organisms:
        find_clusters_in_organism(organism, conserve, gap, require=require)

    return organisms
