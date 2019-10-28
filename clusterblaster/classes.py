#!/usr/bin/env python3

"""
This module stores the classes (Organism, Scaffold, Hit) used in clusterblaster.
"""

import re


def generate_header_string(text, symbol="-"):
    """Generate an underlined header string.
    The underline is the specified symbol repeated for the length of the text.
    """
    length = len(text)
    return f"{text}\n{symbol * length}"


def generate_cluster_table(hits, decimals=4, show_headers=True):
    """Generate human-readable summary of a hit cluster.

    Parameters
    ----------
    hits: list
        Collection of Hit instances.
    decimals: int
        How many decimal points to show.
    show_headers: bool
        Show column headers in output.

    Returns
    -------
    str
        Summary table
    """
    rows = [h.values(decimals) for h in hits]

    headers = [
        "Query",
        "Subject",
        "Identity",
        "Coverage",
        "E-value",
        "Bitscore",
        "Start",
        "End",
        "Strand",
    ]

    if show_headers:
        rows.insert(0, headers)

    # Get lengths of longest values in each column for spacing purposes
    lengths = [max(len(hit[i]) for hit in rows) for i in range(9)]

    # Right-fill each column value with whitespace to width for that column
    return "\n".join(
        "  ".join(f"{hit[i]:{lengths[i]}}" for i in range(9)) for hit in rows
    )


class Organism:
    """The Organism class stores hits on scaffolds"""

    def __init__(self, name, strain):
        self.name = name
        self.strain = strain
        self.scaffolds = {}

    def __str__(self):
        total_scaffolds = len(self.scaffolds)
        total_hits = sum(len(s.hits) for s in self.scaffolds)
        return "{} {} [{} hits on {} scaffolds]".format(
            self.name, self.strain, total_hits, total_scaffolds
        )

    def count_hit_clusters(self):
        """Calculate total amount of hit clusters in this Organism."""
        return sum(len(scaffold.clusters) for scaffold in self.scaffolds.values())

    def summary(self, decimals=4, organism_header=True, scaffold_headers=True):
        """Generate a report of all hit clusters in this Organism."""

        if self.count_hit_clusters() == 0:
            raise ValueError("No hit clusters in this Organism")

        report = "\n\n".join(
            scaffold.summary(decimals=decimals, show_header=scaffold_headers)
            for scaffold in self.scaffolds.values()
            if scaffold.clusters
        )

        if organism_header:
            header = generate_header_string(self.full_name, "=")
            return f"{header}\n{report}"

        return report

    @property
    def full_name(self):
        if self.strain in self.name or not self.strain:
            return f"{self.name}"
        return f"{self.name} {self.strain}"


class Scaffold:
    """The Scaffold object stores Hit objects, as well as the logic for detecting Hit
    clusters.
    """

    def __init__(self, accession, hits=None):
        self.accession = accession
        self.hits = hits if hits else []
        self.clusters = []

    def __str__(self):
        return "{} [{} hits in {} clusters]".format(
            self.accession, len(self.hits), len(self.clusters)
        )

    def summary(self, show_header=True, decimals=4):
        """Generate a summary of hit clusters on this Scaffold."""
        if not self.clusters:
            raise ValueError("No clusters on this Scaffold")

        report = "\n\n".join(
            generate_cluster_table(cluster, decimals=decimals)
            for cluster in self.clusters
        )

        if show_header:
            header = generate_header_string(self.accession)
            return f"{header}\n{report}"

        return report

    @property
    def protein_uids(self):
        return ",".join(hit.target for hit in self.hits)


class Hit:
    """The Hit class stores BLAST hits and their genomic contexts.

    They are first instantiated when parsing BLAST results, and are then updated with
    genomic information after NCBI is queried for IPG.
    """

    def __init__(self, query, subject, identity, coverage, evalue, bitscore):
        self.query = query

        if "gb" in subject or "ref" in subject:
            subject = re.search(r"\|([A-Za-z0-9\._]+)\|", subject).group(1)

        self.subject = subject
        self.bitscore = float(bitscore)
        self.identity = float(identity)
        self.coverage = float(coverage)
        self.evalue = float(evalue)
        self.start = None
        self.end = None

    def __str__(self):
        return "\t".join(self.report())

    def values(self, decimals=4):
        return [
            self.query,
            self.subject,
            f"{self.identity:.{decimals}}",
            f"{self.coverage:.{decimals}}",
            f"{self.evalue:.{decimals}}",
            f"{self.bitscore:.{decimals}}",
            str(self.start),
            str(self.end),
            self.strand,
        ]
