#!/usr/bin/env python3

"""
The `context` module provides methods for interacting with the NCBI's Identical
Protein Group (IPG) resource, grouping of `Hit` objects by genomic scaffold and
organism, and identifying hit clusters.

The main functionality of this module is wrapped up in the `search()` function.

For example, given a list of `Hit` objects, we can query NCBI for genomic coordinates of
each hit, and automatically group them into `Scaffold` and `Organism` objects as follows:


"""


import logging
from collections import defaultdict, namedtuple
from operator import attrgetter

import requests

from cblaster import database
from cblaster.classes import Organism, Scaffold


LOG = logging.getLogger(__name__)


def efetch_IPGs(ids, output_handle=None):
    """Queries the Identical Protein Groups (IPG) resource for given IDs.

    The NCBI caps Efetch requests at 10000 maximum returned records (retmax=10000)
    so this function splits the supplied IDs into chunks of 10000 and queries NCBI
    individually for each chunk.

    Args:
        ids (list): Valid NCBI sequence identifiers.
        output_handle (file handle): File handle to write to.
    Returns:
        List of rows from resulting IPG table, split by newline.
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


def parse_IPG_table(results, hits):
    """Parse IPG results table from `efetch_IPGs`.

    Process:
        1. Parse individual entries in IPG table, grouping by UID
        2. Find representative `Hit` object and make a copy for each IPG entry
        3. Update each `Hit` copy with genomic location from IPG
        4. Save `Hit` copies to corresponding `Organism` objects

    Though this will capture a lot of redundant information (e.g. gene clusters
    that have been deposited separately to their parent genomes) it should not
    miss any hits.

    Args:
        results (list): Results from IPG search.
        hits (list): Hit objects that were used to query NCBI.
    Returns:
        Organism objects containing hits sorted into genomic scaffolds.
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
    for line in results:
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

            # Create new Organism and Scaffold objects on first encounters
            if st in org:
                org = org.replace(st, "").strip()
            if st not in organisms[org]:
                organisms[org][st] = Organism(name=org, strain=st)
            if acc not in organisms[org][st].scaffolds:
                organisms[org][st].scaffolds[acc] = Scaffold(acc)

            # Copy the original Hit object and add contextual information
            h = hit.copy(
                subject=entry.protein_id,
                end=int(entry.end),
                start=int(entry.start),
                strand=entry.strand,
            )
            organisms[org][st].scaffolds[acc].hits.append(h)

    return [organism for strains in organisms.values() for organism in strains.values()]


def find_identifier(qualifiers):
    """Finds an identifier from a dictionary of feature qualifiers.

    This function selects for the following fields in decreasing order:
    protein_id, locus_tag, ID and Gene. This should cover most cases where CDS
    features do not have protein ID's.

    Args:
        qualifiers (dict): Feature qualifiers parsed with genome2json.
    Returns:
        Identifier, if found, otherwise None.
    """
    for field in ("protein_id", "locus_tag", "ID", "Gene"):
        try:
            return qualifiers[field]
        except KeyError:
            pass
    return None


def query_local_DB(hits, database):
    """Build Organisms/Scaffolds using database.DB instance.

    Retrieves corresponding Protein object from the DB by passing it the Hit subject
    attribute, which should be a 4 field '|' delimited header string containing the full
    lineage of the protein (e.g. organism|strain|scaffold|protein).

    Internally uses defaultdict to build up hierarchy of unique organisms, then returns
    list of just each organism dictionary.

    Args:
        hits (list): Hit objects created during cblaster search.
        database (database.DB): cblaster database object.
    Returns:
        Organism objects containing hits sorted into genomic scaffolds.
    """

    organisms = defaultdict(dict)

    for hit in hits:
        # Hit headers should follow form "i_j_k", where i, j and k refer to the
        # database indexes of organisms, scaffolds and proteins, respectively.
        # e.g. >2_56_123 => 123rd protein of 56th scaffold of the 2nd organism
        try:
            i, j, k = [int(index) for index in hit.subject.split("_")]
        except ValueError:
            LOG.exception("Hit has malformed header")

        organism = database.organisms[i]
        scaffold = organism.scaffolds[j]
        protein = scaffold.features[k]

        # For brevity...
        org = organism.name
        st = organism.strain
        sc = scaffold.accession

        # Instantiate new Organism/Scaffold objects on first encounter
        if st not in organisms[org]:
            organisms[org][st] = Organism(org, st)
        if sc not in organisms[org][st].scaffolds:
            organisms[org][st].scaffolds[sc] = Scaffold(sc)

        # Want to report just protein ID, not lineage
        identifier = find_identifier(protein.qualifiers)
        if not identifier:
            LOG.warning("Could not find identifier for hit %, skipping", hit.subject)
            continue
        hit.subject = identifier

        # Save genomic location on the Hit instance
        hit.start = protein.location.min()
        hit.end = protein.location.max()
        hit.strand = protein.location.strand

        organisms[org][st].scaffolds[sc].hits.append(hit)

    return [
        organism
        for strains in organisms.values()
        for organism in strains.values()
    ]


def hits_contain_required_queries(hits, queries):
    """Checks that a group of Hit objects has one Hit per query.

    Args:
        hits (list): Hit objects to be tested.
        queries (list): Names of query sequences to test for.
    Returns:
        True or False
    """
    bools = [False] * len(queries)
    for index, query in enumerate(queries):
        for hit in hits:
            if hit.query == query:
                bools[index] = True
                break
    return all(bools)


def hits_contain_unique_queries(hits, threshold=3):
    """Checks for at least threshold unique queries in a collection of Hit objects."""
    return len(set(hit.query for hit in hits)) >= threshold


def find_clusters(hits, require=None, unique=3, min_hits=3, gap=20000):
    """Finds clusters of Hit objects matching user thresholds.

    Args:
        hits (list): Collection of Hit objects to find clusters in.
        require (list): Names of query sequences that must be represented in a cluster.
        unique (int): Unique query sequence threshold.
        min_hits (int): Minimum number of hits in a hit cluster.
        gap (int): Maximum intergenic distance (bp) between any two hits in a cluster.
    Returns:
        Clusters of Hit objects.
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


def find_clusters_in_organism(organism, unique=3, min_hits=3, gap=20000, require=None):
    """Runs find_clusters() on all scaffolds in an organism."""
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
    """Gets the genomic context for a collection of Hit objects.

    This function wraps all functionality in this module, allowing for Organism objects
    to be directly instantiated from Hit objects insantiated prior.

    >>> organisms = search(
    ...     hits,
    ...     unique=3,  # number of query sequences that MUST be hit in a cluster
    ...     min_hits=3,  # minimum number of hits in a cluster
    ...     gap=20000,  # maximum intergenic gap
    ...     require=["q1", "q2"],  # specific query sequences that MUST be in a cluster
    ...     json_db=None,  # query a JSON database created using cblaster makedb
    ... )

    Args:
        hits (list): Collection of Hit objects to find clusters in.
        require (list): Names of query sequences that must be represented in a cluster.
        unique (int): Unique query sequence threshold.
        min_hits (int): Minimum number of hits in a hit cluster.
        gap (int): Maximum intergenic distance (bp) between any two hits in a cluster.
        json_db (str): Path to a JSON database created with cblaster makedb.
    Returns:
        Dictionary of Organism objects keyed on species name.
    """
    if json_db:
        LOG.info("Loading JSON database: %s", json_db)
        db = database.Database.from_json(json_db)
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
