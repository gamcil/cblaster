#!/usr/bin/env python3

"""
This module stores the classes (Organism, Scaffold, Hit) used in cblaster.
"""

import re


def generate_header_string(text, symbol="-"):
    """Generate an underlined header string.
    The underline is the specified symbol repeated for the length of the text.
    """
    length = len(text)
    return f"{text}\n{symbol * length}"


def generate_cluster_table(hits, decimals=4, show_headers=True, human=True):
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

    if show_headers:
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
        rows.insert(0, headers)

    if human:
        # Get lengths of longest values in each column for spacing purposes
        lengths = [max(len(hit[i]) for hit in rows) for i in range(9)]

        # Right-fill each column value with whitespace to width for that column
        return "\n".join(
            "  ".join(f"{hit[i]:{lengths[i]}}" for i in range(9)) for hit in rows
        )
    return "\n".join(",".join(hit) for hit in rows)


class Organism:
    """The Organism class stores hits on scaffolds"""

    def __init__(self, name, strain):
        self.name = name
        self.strain = strain
        self.scaffolds = {}

    def __str__(self):
        total_scaffolds = len(self.scaffolds)
        total_hits = sum(len(s.hits) for s in self.scaffolds.values())
        return "ORGANISM: {} {} [{} hits on {} scaffolds]".format(
            self.name, self.strain, total_hits, total_scaffolds
        )

    def count_hit_clusters(self):
        """Calculate total amount of hit clusters in this Organism."""
        return sum(len(scaffold.clusters) for scaffold in self.scaffolds.values())

    def summary(self, decimals=4, human=True, headers=True):
        """Generate a report of all hit clusters in this Organism."""

        if self.count_hit_clusters() == 0:
            raise ValueError("No hit clusters in this Organism")

        report = "\n\n".join(
            scaffold.summary(decimals=decimals, human=human, show_header=headers)
            for scaffold in self.scaffolds.values()
            if scaffold.clusters
        )

        if headers:
            header = generate_header_string(self.full_name, "=")
            return f"{header}\n{report}"

        return report

    @property
    def full_name(self):
        if not self.strain or self.strain in self.name:
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
        return "SCAFFOLD: {} [{} hits in {} clusters]".format(
            self.accession, len(self.hits), len(self.clusters)
        )

    def summary(self, human=True, show_header=True, decimals=4):
        """Generate a summary of hit clusters on this Scaffold."""
        if not self.clusters:
            raise ValueError("No clusters on this Scaffold")

        report = "\n\n".join(
            generate_cluster_table(
                cluster, decimals=decimals, show_headers=show_header, human=human
            )
            for cluster in self.clusters
        )

        if show_header:
            header = generate_header_string(self.accession)
            return f"{header}\n{report}"

        return report


class Hit:
    """The Hit class stores BLAST hits and their genomic contexts.

    They are first instantiated when parsing BLAST results, and are then updated with
    genomic information after NCBI is queried for IPG.
    """

    def __init__(
        self,
        query,
        subject,
        identity,
        coverage,
        evalue,
        bitscore,
        start=None,
        end=None,
        strand=None,
    ):
        self.query = query

        if "gb" in subject or "ref" in subject:
            subject = re.search(r"\|([A-Za-z0-9\._]+)\|", subject).group(1)

        self.subject = subject
        self.bitscore = float(bitscore)
        self.identity = float(identity) / 100
        self.coverage = float(coverage) / 100
        self.evalue = float(evalue)

        self.start = int(start) if start is not None else None
        self.end = int(end) if end is not None else None
        self.strand = strand

    def __str__(self):
        return (
            f"Hit: {self.query} - {self.subject}:"
            f"{self.start}-{self.end}[{self.strand}] "
            f" {self.identity:.2%}/{self.coverage:.2%}"
        )

    def values(self, decimals=4):
        """Format all attributes of this hit for printing.

        Parameters
        ----------
        decimals: int
            Maximum number of decimal points to report. Note that this applies to
            rounding, not necessarily display; i.e. if the attribute is 50 and
            decimals=4, this function will return '50', not '50.0000'. The e-value is
            also capped using this parameter, but using exponent formatting.

        Returns
        -------
        list
            All attributes of this Hit, formatted as str.
        """
        return [
            self.query,
            self.subject,
            f"{round(self.identity, decimals):g}",
            f"{round(self.coverage, decimals):g}",
            f"{self.evalue:.{decimals}g}",
            f"{round(self.bitscore, decimals):g}",
            str(self.start),
            str(self.end),
            self.strand,
        ]
