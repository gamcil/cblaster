#!/usr/bin/env python3

"""
This module provides methods for retrieving contextual information about search hits
from the NCBI's Identical Protein Group (IPG) resource and JSON databases created using
cblaster makedb. It also provides methods for filtering and clustering hits based on
user specified thresholds.

The main functionality of this module is wrapped up in the search() function.
Given a list of Hit objects (resulting from a cblaster search), we can call it:

>>> search(hits)

This will:
1. Search the identifiers of the subject hits against the IPG and retrieve results
2. Parse protein groups from the IPG table
3. Create Subject objects for each entry in any given IPG, containing copies of Hit objects, grouped by organism and scaffold
4. Identify clusters based on user thresholds for intergenic distance, copy number, etc
5. De-duplicate clusters within an organism, where all Subject objects in any two clusters are members of the same IPG

Note that cblaster uses Subject objects, not Hit objects, for this step. A Subject
refers to a unique protein in any given organism, which can be hit in a search any
number of times. A Hit object just stores information about one of these hits (query
and subject protein identifiers, score values). In order to prevent redundancy, Subjects
are created for each unique protein found in the search, which then store the Hit
objects linking them to the query sequences.
"""


import logging
from collections import defaultdict, namedtuple
from itertools import combinations, product
from operator import attrgetter
from functools import partial

import requests
import numpy as np

from cblaster import database
from cblaster.classes import Organism, Scaffold, Subject


LOG = logging.getLogger(__name__)


def efetch_IPGs(ids, output_file=None):
    """Queries the Identical Protein Groups (IPG) resource for given IDs.

    The NCBI caps Efetch requests at 10000 maximum returned records (retmax=10000)
    so this function splits the supplied IDs into chunks of 10000 and queries NCBI
    individually for each chunk.

    Args:
        ids (list): Valid NCBI sequence identifiers.
        output_file (str): File to write the ipg table to.
    Returns:
        List of rows from resulting IPG table, split by newline.
    """

    # Split into chunks since retmax=10000
    chunks = [ids[i: i + 10000] for i in range(0, len(ids), 10000)]

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
                f"Error fetching sequences from NCBI [code {response.status_code}]."
            )

        table += response.text

    if output_file:
        LOG.info("Writing IPG table to %s", output_file)
        with open(output_file, "w") as f:
            f.write(table)

    return table.split("\n")


def parse_IP_groups(results):
    """Parse groups from an Identical Protein Groups (IPG) table.

    This function converts rows in the IPG table to namedtuple objects which have
    attributes corresponding to each field in the table. These objects are grouped
    by the IPG they belong to.

    Args:
        results (list): Rows in the IPG table.
    Returns:
        Dictionary of table entries (namedtuple objects) grouped by IPG.
    """
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
        if not line or line.startswith("Id\tSource") or line.isspace() \
                or "skipping" in line:
            continue
        ipg, *fields = line.strip("\n").split("\t")
        entry = Entry(*fields)
        groups[ipg].append(entry)
    return groups


def find_IPG_hits(group, hit_dict):
    """Finds all hits to query sequences in an identical protein group.

    Args:
        group (list): Entry namedtuples from an IPG from parse_IP_groups().
        hit_dict (dict): All Hit objects found during cblaster search.
    Returns:
        List of all Hit objects corresponding to any proteins in a IPG.
    """
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


def group_hits(hits):
    hit_dict = defaultdict(list)
    for hit in hits:
        hit_dict[hit.subject].append(hit)
    return hit_dict


def parse_IPG_table(results, hits):
    """Links Hit objects to their genomic context from an IPG table.

    This function:
    1. Parses entries from the table, grouped by IPG number with parse_IP_groups()
    2. For each group:

    a) Find Hit objects linked to any member of the group with find_IPG_hits()
    b) For each group member, create a Subject object, then place it on its
       corresponding Scaffold and Organism objects (creating new objects when new
       scaffolds and organisms are encountered)
    c) Add Hit objects to every Subject object in the group

    Args:
        results (list): Results from IPG search.
        hits (list): Hit objects that were used to query NCBI.
    Returns:
        Organism objects containing hits sorted into genomic scaffolds.
    """

    # Group hits by their subject IDs
    hit_dict = group_hits(hits)

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
            # Avoid incomplete entries
            if not all([entry.scaffold, entry.start, entry.end, entry.strand]):
                continue

            # Avoid vectors, single Gene nucleotide entries, etc
            if entry.source not in ("RefSeq", "INSDC",):
                continue

            # Test unique scaffold and coordinates - sometimes will have many identical
            # protein_id in one CDS region
            test = (entry.scaffold, entry.start, entry.end)
            if test in seen:
                continue
            seen.add(test)

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
                hits=[hit.copy(subject=entry.protein_id) for hit in hit_list],
                name=entry.protein_id,
                ipg=ipg,
                end=int(entry.end),
                start=int(entry.start),
                strand=entry.strand,
            )

            organisms[org][st].scaffolds[acc].subjects.append(subject)

    return [
        organism
        for strains in organisms.values()
        for organism in strains.values()
    ]


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


def query_local_DB(hits, db):
    """Queries a local SQLite3 database created using the makedb module.
    """
    organisms = {}
    hit_dict = defaultdict(list)
    for hit in hits:
        hit_dict[hit.subject].append(hit)
    for (
        rowid,
        name,
        start_pos,
        end_pos,
        strand,
        scaffold,
        organism
    ) in database.query_genes(list(hit_dict), db):
        if organism not in organisms:
            organisms[organism] = Organism(organism, "")
        if scaffold not in organisms[organism].scaffolds:
            organisms[organism].scaffolds[scaffold] = Scaffold(scaffold)
        hits = hit_dict[str(rowid)]
        for hit in hits:
            hit.subject = name
        subject = Subject(
            id=rowid,
            name=name,
            hits=hits,
            start=int(start_pos),
            end=int(end_pos),
            strand="+" if strand == 1 else "-"
        )
        organisms[organism].scaffolds[scaffold].subjects.append(subject)
    return [organism for organism in organisms.values()]


def cluster_satisfies_conditions(cluster, require=None, unique=3, minimum=3):
    """Tests if a cluster of Subjects meets query conditions.

    Finds all unique query sequences in hits, then returns True if total number is above
    unique threshold, and any required queries are represented.
    """
    queries = set(hit.query for subject in cluster for hit in subject.hits)
    return (
        len(cluster) >= minimum
        and len(queries) >= unique
        and (queries.issuperset(require) if require else True)
    )


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
            return [subjects]
        return []

    sorted_subjects = sorted(subjects, key=attrgetter("start"))
    first = sorted_subjects.pop(0)
    group, border = [first], first.end

    rules_satisfied = partial(
        cluster_satisfies_conditions,
        require=require,
        unique=unique,
        minimum=min_hits
    )

    for subject in sorted_subjects:
        if subject.start <= border + gap:
            group.append(subject)
            border = max(border, subject.end)
        else:
            if rules_satisfied(group):
                yield group
            group, border = [subject], subject.end
    if rules_satisfied(group):
        yield group


def clusters_are_identical(one, two):
    """Tests if two collections of Subject objects are identical.

    This function compares clusters element-wise for IPG numbers, returning True
    if all are identical.
    """
    if not len(one) == len(two):
        return False
    for subA, subB in zip(one, two):
        if not subA.ipg or not subB.ipg:
            return False
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
            if clusters_are_identical(one, two):
                remove[scafB.accession].append(two)
    for accession, clusters in remove.items():
        scaffold = organism.scaffolds[accession]
        scaffold.clusters = [c for c in scaffold.clusters if c not in clusters]


def find_clusters_in_organism(
    organism,
    unique=3,
    min_hits=3,
    gap=20000,
    require=None,
    remote=True,
    query_sequence_order=None
):
    """Runs find_clusters() on all scaffolds in an organism."""
    for scaffold in organism.scaffolds.values():
        clusters = find_clusters(
            scaffold.subjects,
            unique=unique,
            min_hits=min_hits,
            gap=gap,
            require=require,
        )
        scaffold.add_clusters(clusters, query_sequence_order=query_sequence_order)
        LOG.debug(
            "Organism: %s, Scaffold: %s, Clusters: %i",
            organism.full_name,
            scaffold.accession,
            len(scaffold.clusters),
        )
    if remote:
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
                if len(subject.hits) == 0:
                    scaffold.remove_subject(subject)
            clusters = find_clusters(
                scaffold.subjects,
                gap=gap,
                min_hits=min_hits,
                require=require,
                unique=unique,
            )
            scaffold.clusters = []

            scaffold.add_clusters(clusters, query_sequence_order=session.queries)
        deduplicate(organism)


def calculate_gne(session):
    """
    No. clusters
    Mean gn size
    Median gn size
    """
    clusters = [
        cluster.end - cluster.start
        for organism in session.organisms
        for accession, scaffold in organism.scaffolds.items()
        for cluster in scaffold.clusters
    ]
    if not clusters:
        return 0, 0, 0
    return (
        len(clusters),
        float(np.mean(clusters)),
        int(np.median(clusters))
    )


def estimate_neighbourhood(session, max_gap=100000, samples=100, scale="linear"):
    """Estimate gene neighbourhood of a cblaster session."""
    if scale == "linear":
        space = np.linspace(0, max_gap, num=samples)
    elif scale == "log":
        space = np.geomspace(1, max_gap, num=samples)
    else:
        raise ValueError("Invalid scale specified, expected 'linear' or 'log'")
    results = []
    for value in space:
        value = int(value)
        filter_session(session, gap=value)
        clusters, means, medians = calculate_gne(session)
        result = {
            "gap": value,
            "means": means,
            "medians": medians,
            "clusters": clusters,
        }
        results.append(result)
    return results


def search(
    hits,
    sqlite_db=None,
    unique=3,
    min_hits=3,
    gap=20000,
    require=None,
    ipg_file=None,
    query_sequence_order=None,
):
    """Gets the genomic context for a collection of Hit objects.

    This function wraps all functionality in this module, allowing for Organism objects
    to be directly instantiated from Hit objects insantiated prior.

    >>> organisms = search(
    ...     hits,
    ...     unique=3,  # number of query sequences that MUST be hit in a cluster
    ...     min_hits=3,  # minimum number of hits in a cluster
    ...     gap=20000,  # maximum intergenic gap
    ...     require=["q1", "q2"],  # specific query sequences that MUST be in a cluster
    ... )

    Args:
        hits (list): Collection of Hit objects to find clusters in.
        sqlite_db (str): path to sqlite database.
        require (list): Names of query sequences that must be represented in a cluster.
        unique (int): Unique query sequence threshold.
        min_hits (int): Minimum number of hits in a hit cluster.
        gap (int): Maximum intergenic distance (bp) between any two hits in a cluster.
        query_sequence_order (list): list of sequences of the order in the query file, is
        ipg_file (str): file to save the ipg table into.
    Returns:
        Dictionary of Organism objects keyed on species name.
    """
    if sqlite_db:
        LOG.info("Querying local SQLite3 database: %s", sqlite_db)
        organisms = query_local_DB(hits, sqlite_db)
    else:
        rows = efetch_IPGs(
            [hit.subject for hit in hits],
            output_file=ipg_file
        )
        organisms = parse_IPG_table(rows, hits)

    LOG.info("Searching for clustered hits across %i organisms", len(organisms))
    for organism in organisms:
        find_clusters_in_organism(
            organism,
            unique=unique,
            min_hits=min_hits,
            gap=gap,
            require=require,
            query_sequence_order=query_sequence_order,
            remote=sqlite_db is None,
        )

    return organisms
