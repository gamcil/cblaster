#!/usr/bin/env python3

"""
This module stores the classes (Organism, Scaffold, Hit) used in cblaster.
"""

import re
import json


def generate_header_string(text, symbol="-"):
    """Generates a 2-line header string with underlined text.

    >>> header = generate_header_string("header string", symbol="*")
    >>> print(header)
    header string
    *************
    """
    return f"{text}\n{symbol * len(text)}"


def generate_cluster_table(hits, decimals=4, show_headers=True, human=True):
    """Generates a summary table for a hit cluster.

    Args:
        hits (list): collection of Hit objects
        decimals (int): number of decimal points to show
        show_headers (bool): show column headers in output
        human (bool): use human-readable format
    Returns:
        summary table
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


def count_query_hits(queries, hits):
    """Counts total hits per query in a colllection of `Hit` objects.

    >>> queries = ["query1", "query2", "query3"]
    >>> hits = [
    ...     Hit(query="query1"),
    ...     Hit(query="query2"),
    ...     Hit(query="some_other_query"),
    ... ]
    >>> count_query_hits(queries, hits)
    [1, 1, 0]

    Args:
        hits (list): Hit objects
    Returns:
        List of per-query counts corresponding to input.
    """
    return [sum(query == hit.query for hit in hits) for query in queries]


def get_max_hit_identities(queries, hits):
    """Get the maximum hit identity per query in a collection of `Hit` objects.

    >>> queries = ["query1", "query2", "query3"]
    >>> hits = [
    ...     Hit(query="query1", identity=0.9),
    ...     Hit(query="query2", identity=0.4),
    ...     Hit(query="query2", identity=0.7),
    ... ]
    >>> get_max_hit_identities(queries, hits)
    [0.9, 0.7, 0]
    """
    return [
        max([hit.identity if query == hit.query else 0 for hit in hits])
        for query in queries
    ]


class Serializer:
    """JSON serialisation mixin class.

    Classes that inherit from this class should implement `to_dict` and
    `from_dict` methods.
    """

    def to_dict(self):
        """Serialises class to dict."""
        raise NotImplementedError

    @classmethod
    def from_dict(self, d):
        """Loads class from dict."""
        raise NotImplementedError

    def to_json(self, fp=None, **kwargs):
        """Serialises class to JSON."""
        d = self.to_dict()
        if fp:
            json.dump(d, fp, **kwargs)
        else:
            return json.dumps(d, **kwargs)

    @classmethod
    def from_json(cls, js):
        """Instantiates class from JSON handle."""
        if isinstance(js, str):
            d = json.loads(js)
        else:
            d = json.load(js)
        return cls.from_dict(d)


class Session(Serializer):
    """Stores the state of a cblaster search.

    This class stores query proteins, search parameters, Organism objects created
    during searches, as well as methods for generating summary tables. It can also be
    dumped to/loaded from JSON for re-filtering, plotting, etc.

    >>> s = Session()
    >>> with open("session.json", "w") as fp:
    ...     s.to_json(fp)
    >>> with open("session.json") as fp:
    ...     s2 = Session.from_json(fp)
    >>> s == s2
    True

    Attributes:
        queries (list): Names of query sequences.
        params (dict): Search parameters.
        organisms (list): Organism objects created in a search.
    """

    def __init__(self, queries, params, organisms=None):
        self.queries = queries
        self.params = params
        self.organisms = organisms if organisms else []

    def to_dict(self):
        return {
            "queries": self.queries,
            "params": self.params,
            "organisms": [o.to_dict() for o in self.organisms],
        }

    @classmethod
    def from_dict(cls, d):
        return cls(
            d["queries"],
            d["params"],
            organisms=[Organism.from_dict(o) for o in d["organisms"]],
        )

    def form_matrices(self, html=False):
        """Form 2D count matrices required for plotting.

        If `html` is False, generated dendrogram leaves will contain LaTeX formatting.
        """

        def form_row(organism, accession, cluster):
            try:
                genus, species = organism.name.split(" ", 1)
                if html:
                    name = f"{genus} {species} {organism.strain}"
                else:  # matplotlib with LaTeX formatting
                    name = f"$\it{{{genus[0]}. {species}}}$ {organism.strain}"
            except ValueError:
                name = organism.name
            scaf = f"{accession}:{cluster[0].start}-{cluster[-1].end}"
            cnts = self.count_query_hits(cluster)
            idts = self.get_max_hit_identities(cluster)
            return name, scaf, cnts, idts

        names = []  # Species names for dendrogram yticklabels
        scafs = []  # Scaffold positions for matshow yticklabels
        counts = []  # Total # of hits in cluster per query
        idents = []  # Maximum identity of single hit per query

        for organism in self.organisms:
            for accession, scaffold in organism.scaffolds.items():
                for cluster in scaffold.clusters:
                    name, scaf, cnts, idts = form_row(organism, accession, cluster)
                    names.append(name)
                    scafs.append(scaf)
                    counts.append(cnts)
                    idents.append(idts)

        return names, scafs, counts, idents

    def _summary(self, human=True, headers=True):
        """Generate summary of >1 Organisms, print to console or write to file.

        Args:
            human (bool): Use human-readable format.
            headers (bool): Show table headers.
        Returns:
            The summary table.
        """
        return "\n\n\n".join(
            organism.summary(headers=headers, human=human)
            for organism in self.organisms
            if organism.total_hit_clusters > 0
        )

    def _binary(self, human=False, headers=True, identity=False):
        """Generate a binary summary table.

        For example:

        Organism  Scaffold  Start  End    Query1  Query2  Query3  Query4
        Org 1     Scaf_1    1      20000  2       1       1       1
        Org 1     Scaf_3    3123   40302  0       1       0       1
        """

        columns = len(self.queries) + 4

        rows = [
            [
                organism.full_name,
                accession,
                str(cluster[0].start),
                str(cluster[-1].end),
                *[
                    str(x)
                    for x in (
                        get_max_hit_identities(self.queries, cluster)
                        if identity
                        else count_query_hits(self.queries, cluster)
                    )
                ],
            ]
            for organism in self.organisms
            for accession, scaffold in organism.scaffolds.items()
            for cluster in scaffold.clusters
        ]

        if headers:
            rows.insert(0, ["Organism", "Scaffold", "Start", "End", *self.queries])

        if human:
            # Calculate lengths of each column for spacing
            lengths = [max(len(row[i]) for row in rows) for i in range(columns)]
            table = "\n".join(
                "  ".join(f"{row[i]:{lengths[i]}}" for i in range(columns))
                for row in rows
            )
        else:
            table = "\n".join(",".join(row) for row in rows)

        return table

    def format(self, form, fp, human=True, headers=True, **kwargs):
        """Generates a summary table.

        Args:
            form (str): Type of table to generate ('summary' or 'binary').
            fp (file handle): File handle to write to.
            human (bool): Use human-readable format.
            headers (bool): Show table headers.
        Raises:
            ValueError: `form` not 'binary' or 'summary'
        Returns:
            Summary table.
        """
        if form == "summary":
            table = self._summary(human=human, headers=headers, **kwargs)
        elif form == "binary":
            table = self._binary(human=human, headers=headers, **kwargs)
        else:
            raise ValueError("Expected 'summary' or 'binary'")
        print(table, file=fp)


class Organism(Serializer):
    """A unique organism containing hits found in a cblaster search.

    Every strain (or lack thereof) is a unique Organism, and will be reported
    separately in cblaster results.

    Attributes:
        name (str): Organism name, typically the genus and species epithet.
        strain (str): Strain name of this organism, e.g. CBS 536.65.
        scaffolds (dict): Scaffold objects belonging to this organism.
    """

    def __init__(self, name, strain, scaffolds=None):
        self.name = name
        self.strain = strain
        self.scaffolds = scaffolds if scaffolds else {}

    def __str__(self):
        total_scaffolds = len(self.scaffolds)
        total_hits = sum(len(s.hits) for s in self.scaffolds.values())
        return "ORGANISM: {} {} [{} hits on {} scaffolds]".format(
            self.name, self.strain, total_hits, total_scaffolds
        )

    @property
    def total_hit_clusters(self):
        """Counts total amount of hit clusters in this Organism."""
        return sum(len(scaffold.clusters) for scaffold in self.scaffolds.values())

    def summary(self, decimals=4, human=True, headers=True):
        """Generates a summary table of the organism.

        Args:
            decimals (int): Total decimal places to show in score values.
            human (bool): Use human-readable format.
            headers (bool): Show table headers.
        Returns:
            The summary table.
        """

        if self.total_hit_clusters == 0:
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
        """The full name (including strain) of the organism.
        Note: if strain found in name, returns just name.
        """
        if not self.strain or self.strain in self.name:
            return f"{self.name}"
        return f"{self.name} {self.strain}"

    def to_dict(self):
        return {
            "name": self.name,
            "strain": self.strain,
            "scaffolds": [scaffold.to_dict() for scaffold in self.scaffolds.values()],
        }

    @classmethod
    def from_dict(cls, d):
        return cls(
            d["name"],
            d["strain"],
            scaffolds={
                scaffold["accession"]: Scaffold.from_dict(scaffold)
                for scaffold in d["scaffolds"]
            },
        )


class Scaffold(Serializer):
    """A genomic scaffold containing hits found in a cblaster search.

    Attributes:
        accession (str): Name of this scaffold, typically NCBI accession.
        hits (list): Hit objects located on this scaffold.
        clusters (list): Clusters of hits identified on this scaffold.
    """

    def __init__(self, accession, clusters=None, hits=None):
        self.accession = accession
        self.hits = hits if hits else []
        self.clusters = clusters if clusters else []

    def __str__(self):
        return "SCAFFOLD: {} [{} hits in {} clusters]".format(
            self.accession, len(self.hits), len(self.clusters)
        )

    def summary(self, human=True, show_header=True, decimals=4):
        """Generates a summary of hit clusters on this Scaffold.

        Args:
            human (bool): Use human-readable format.
            show_header (bool): Show table headers.
            decimals (int): Total decimal places to show in score values.
        Returns:
            The summary table.
        """
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

    def to_dict(self):
        return {
            "accession": self.accession,
            "hits": [hit.to_dict() for hit in self.hits],
            "clusters": [
                [self.hits.index(hit) for hit in cluster] for cluster in self.clusters
            ],
        }

    @classmethod
    def from_dict(cls, d):
        hits = [Hit.from_dict(hit) for hit in d["hits"]]
        clusters = [[hits[ix] for ix in cluster] for cluster in d["clusters"]]
        return cls(accession=d["accession"], hits=hits, clusters=clusters)


class Hit(Serializer):
    """A BLAST hit identified during a cblaster search.

    This class is first instantiated when parsing BLAST results, and is then updated
    with genomic coordinates after querying either the Identical Protein Groups (IPG)
    resource on NCBI, or a local JSON database.

    Attributes:
        query (str): Name of query sequence.
        subject (str): Name of subject sequence.
        identity (float): Percentage identity (%) of hit.
        coverage (float): Query coverage (%) of hit.
        evalue (float): E-value of hit.
        bitscore (float): Bitscore of hit.
        start (int): Start of subject sequence on corresponding scaffold.
        end (int): End of subject sequence on corresponding scaffold
        strand (str): Orientation of subject sequence ('+' or '-').
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
        self.identity = float(identity)
        self.coverage = float(coverage)
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

    def copy(self, **kwargs):
        """Creates a copy of this Hit with any additional args."""
        copy = Hit(**self.__dict__)
        for key, val in kwargs.items():
            setattr(copy, key, val)
        return copy

    def values(self, decimals=4):
        """Formats hit attributes for printing.

        Args:
            decimals (int): Total decimal places to show in score values.
        Returns:
            List of formatted attribute strings.
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

    def to_dict(self):
        return {
            "query": self.query,
            "subject": self.subject,
            "identity": self.identity,
            "coverage": self.coverage,
            "evalue": self.evalue,
            "bitscore": self.bitscore,
            "start": self.start,
            "end": self.end,
            "strand": self.strand,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(**d)
