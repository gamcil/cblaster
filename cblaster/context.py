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
from itertools import combinations, product
from operator import attrgetter

import requests

from cblaster import database
from cblaster.classes import Organism, Scaffold, Subject


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


def parse_IP_groups(results):
    """Parse groups from an Identical Protein Groups (IPG) table."""
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
    groups = defaultdict(list)
    for line in results:
        if not line or line.startswith("Id\tSource") or line.isspace():
            continue
        ipg, *fields = line.split("\t")
        entry = Entry(*fields)

        # Avoid incomplete entries
        if not all([entry.scaffold, entry.start, entry.end, entry.strand]):
            continue

        # Avoid vectors, single Gene nucleotide entries, etc
        if entry.source not in ("RefSeq", "INSDC",):
            continue

        groups[ipg].append(entry)
    return groups


def find_IPG_hits(group, hit_dict):
    """Finds all hits to query sequences in an identical protein group."""
    seen = set()
    hits = []
    for entry in group:
        try:
            new = hit_dict.pop(entry.protein_id)
        except KeyError:
            continue
        for h in new:
            if h.query in seen:
                continue
            seen.add(h.query)
            hits.append(h)
    return hits


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

    # Group hits by their subject ID's
    hit_dict = defaultdict(list)
    for hit in hits:
        hit_dict[hit.subject].append(hit)

    # Parse IPGs from the table
    groups = parse_IP_groups(results)

    seen = set()
    organisms = defaultdict(dict)
    for ipg in list(groups):
        group = groups.pop(ipg)

        # Find any hits corresponding to this IPG
        hit_list = find_IPG_hits(group, hit_dict)

        if not hit_list:
            LOG.warning("Found no hits for IPG %s", ipg)
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
            subject = Subject(
                hits=hit_list,
                ipg=group,
                end=int(entry.end),
                start=int(entry.start),
                strand=entry.strand,
            )
            organisms[org][st].scaffolds[acc].subjects.append(subject)

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

    # Form non-redundant dictionary of hits. Each key will become a unique Subject.
    hit_dict = defaultdict(list)
    for hit in hits:
        hit_dict[hit.subject].append(hit)

    for hit_index, hits in hit_dict.items():
        # Hit headers should follow form "i_j_k", where i, j and k refer to the
        # database indexes of organisms, scaffolds and proteins, respectively.
        # e.g. >2_56_123 => 123rd protein of 56th scaffold of the 2nd organism
        try:
            i, j, k = [int(index) for index in hit_index.split("_")]
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
        for hit in hits:
            hit.subject = identifier

        # Save genomic location on the Hit instance
        subject = Subject(
            hits=hits,
            start=protein.location.min(),
            end=protein.location.max(),
            strand=protein.location.strand
        )

        organisms[org][st].scaffolds[sc].hits.append(hit)

    return [
        organism
        for strains in organisms.values()
        for organism in strains.values()
    ]


def cluster_satisfies_conditions(cluster, require=None, unique=3, minimum=3):
    """Tests if a cluster of Subjects meets query conditions.

    Finds all unique query sequences in hits, then returns True if total number is above
    unique threshold, and any required queries are represented.
    """
    required = [False] * (len(require) if require else 0)
    queries = set(hit.query for subject in cluster for hit in subject.hits)
    return len(queries) >= unique and queries.issuperset(required)


def find_clusters(subjects, require=None, unique=3, min_hits=3, gap=20000):
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

    total_subjects = len(subjects)

    if total_subjects < unique:
        return []

    if total_subjects == 1:
        if unique == 1 or min_hits == 1:
            return [hits]
        return []

    sorted_subjects = sorted(subjects, key=attrgetter("start"))
    first = sorted_subjects.pop(0)
    group, border = [first], first.end

    for subject in sorted_subjects:
        if subject.start <= border + gap:
            group.append(subject)
            border = max(border, subject.end)
        else:
            if cluster_satisfies_conditions(group):
                yield group
            group, border = [subject], subject.end

    if cluster_satisfies_conditions(group):
        yield group


def clusters_are_identical(one, two):
    if not len(one) == len(two):
        return False
    for subA, subB in zip(one, two):
        if subA.ipg != subB.ipg:
            return False
    return True


def deduplicate(organism):
    """Removes any duplicate clusters within an Organism.

    Some redundancy is unavoidable due to duplicate entries on NCBI. This function
    attempts to remedy this partially by searching for identical clusters within a
    single, unique Organism. Clusters are compared from start to finish, and are tagged
    for removal if every Subject is of the same IPG for the length of the clusters.
    """
    remove = defaultdict(list)

    for scafA, scafB in combinations(organism.scaffolds.values(), 2):
        for one, two in product(scafA.clusters, scafB.clusters):
            if one == two:
                continue
            if clusters_are_identical(one, two):
                remove[scafB.accession].append(two)

    for accession, clusters in remove.items():
        scaffold = organism.scaffolds[accession]
        scaffold.clusters = [c for c in scaffold.clusters if c not in clusters]


def find_clusters_in_organism(organism, unique=3, min_hits=3, gap=20000, require=None):
    """Runs find_clusters() on all scaffolds in an organism."""
    for scaffold in organism.scaffolds.values():
        scaffold.clusters = list(
            find_clusters(
                scaffold.subjects,
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

    deduplicate(organism)



def filter_session(
    session,
    min_identity=30,
    min_coverage=50,
    max_evalue=0.01,
    gap=20000,
    unique=3,
    min_hits=3,
    require=None,
):
    """Filter a Session object with new thresholds.

    This function is destructive!
    """
    for organism in session.organisms:
        for scaffold in organism.scaffolds.values():
            for subject in scaffold.subjects:
                subject.hits = [
                    hit
                    for hit in subject.hits
                    if (
                        hit.identity > min_identity
                        and hit.coverage > min_coverage
                        and hit.evalue < max_evalue
                    )
                ]
            scaffold.clusters = list(
                find_clusters(
                    scaffold.subjects,
                    gap=gap,
                    min_hits=min_hits,
                    require=require,
                    unique=unique,
                )
            )
        context.deduplicate(organism)



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
