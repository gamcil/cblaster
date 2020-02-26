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
    """Query Identical Protein Groups (IPG) with query sequence IDs.

    The NCBI caps Efetch requests at 10000 maximum returned records (i.e. retmax=10000).
    Thus, this function splits `ids` into chunks of 10000 and queries NCBI individually
    for each chunk.

    Parameters
    ----------
    ids: list, tuple
        NCBI sequence identifiers
    output_handle: open file handle
        File handle to write IPG table to

    Returns
    -------
    table: list
        Results from IPG search split by newline
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
    """Parse IPG results table from `efetch_IPGs`.

    This function first parses individual entries in the IPG table and groups them by
    IPG uid. For each IPG, a representative Hit object is found in `hits` and a copy
    is made for each protein in the IPG. Each copy is updated with the genomic origin
    listed in the IPG table, and then added to the corresponding Organism object (or a
    new Organism if one isn't found). `cblaster` should therefore capture every possible
    cluster on NCBI.

    Parameters
    ----------
    results_handle: open file handle
        Results from identical protein groups search
    hits: list
        Hit objects that were used to query NCBI

    Returns
    -------
    organisms: list
        Organism objects
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

            if st in org:
                org = org.replace(st, "").strip()

            if st not in organisms[org]:
                organisms[org][st] = Organism(name=org, strain=st)

            if acc not in organisms[org][st].scaffolds:
                organisms[org][st].scaffolds[acc] = Scaffold(acc)

            h = hit.copy(
                subject=entry.protein_id,
                end=int(entry.end),
                start=int(entry.start),
                strand=entry.strand,
            )

            organisms[org][st].scaffolds[acc].hits.append(h)

    return [organism for strains in organisms.values() for organism in strains.values()]


def query_local_DB(hits, database):
    """Build Organisms/Scaffolds using database.DB instance.

    Retrieves corresponding Protein object from the DB by passing it the Hit subject
    attribute, which should be a 4 field '|' delimited header string containing the full
    lineage of the protein (e.g. organism|strain|scaffold|protein).

    Internally uses defaultdict to build up hierarchy of unique organisms, then returns
    list of just each organism dictionary.

    Parameters
    ----------
    hits: list
        Hit objects created during cblaster search
    database: database.DB
        cblaster database object

    Returns
    -------
    list
        Organism objects
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


def hits_contain_unique_queries(hits, threshold=3):
    """Check that a group of Hit objects belong to some threshold of unique queries."""
    return len(set(hit.query for hit in hits)) >= threshold


def find_clusters(hits, require=None, unique=3, min_hits=3, gap=20000):
    """Find clustered Hit objects.

    Parameters
    ----------
    hits: list
        Hit objects
    require: iterable
        Names of query sequences that must be represented in a hit cluster
    unique: int
        Minimum number of unique queries represented in a hit cluster
    min_hits: int
        Minimum number of hits in a cluster
    gap: int
        Maximum intergenic distance (bp) between any two hits in a cluster

    Returns
    -------
    groups: list
        Groups of Hit objects.
    """
    if unique < 0 or min_hits < 0 or gap < 0:
        raise ValueError("Expected positive integer")

    total_hits = len(hits)

    if total_hits < unique:
        return []

    if total_hits == 1:
        if unique == 1 or min_hits == 1:
            return [hits]
        return []

    def conditions_met(group):
        req = hits_contain_required_queries(group, require) if require else True
        con = hits_contain_unique_queries(group, unique)
        siz = len(group) >= min_hits
        return req and con and siz

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


def find_clusters_in_organism(
    organism, unique=3, min_hits=3, gap=20000, require=None
):
    """Run find_clusters() on all Hits on Scaffolds in an Organism instance."""
    for scaffold in organism.scaffolds.values():
        scaffold.clusters = list(
            find_clusters(
                scaffold.hits,
                unique=unique,
                min_hits=min_hits,
                gap=gap,
                require=require,
            )
        )
        LOG.debug(
            "Organism: %s, Scaffold: %s, Clusters: %i",
            organism.full_name,
            scaffold.accession,
            len(scaffold.clusters),
        )


def search(hits, unique=3, min_hits=3, gap=20000, require=None, json_db=None):
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
        find_clusters_in_organism(
            organism, unique=unique, min_hits=min_hits, gap=gap, require=require
        )

    return organisms
