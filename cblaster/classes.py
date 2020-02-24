#!/usr/bin/env python3

"""
This module stores the classes (Organism, Scaffold, Hit) used in cblaster.
"""

import re
import json


def generate_header_string(text, symbol="-"):
    """Generate an underlined header string.

    The underline is the specified symbol repeated for the length of the text.
    """
    length = len(text)
    return f"{text}\n{symbol * length}"


def generate_cluster_table(hits, decimals=4, show_headers=True, human=True):
    """Generate a summary table for a hit cluster.

    Parameters
    ----------
    hits: list
        Collection of Hit instances
    decimals: int
        How many decimal points to show
    show_headers: bool
        Show column headers in output
    human: bool
        Generate human readable, not delimited, table

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


class Serializer:
    """JSON serialisation mixin class.

    Classes that inherit from this class should implement `to_dict` and `from_dict`
    methods.
    """

    def to_dict(self):
        raise NotImplementedError

    @classmethod
    def from_dict(self, d):
        raise NotImplementedError

    def to_json(self, fp=None, **kwargs):
        d = self.to_dict()
        if fp:
            json.dump(d, fp, **kwargs)
        else:
            return json.dumps(d, **kwargs)

    @classmethod
    def from_json(cls, js):
        if isinstance(js, str):
            d = json.loads(js)
        else:
            d = json.load(js)
        return cls.from_dict(d)


class Session(Serializer):
    """A search session.

    This class is used to manage state during a cblaster search. It stores IDs of query
    proteins, parameters used, RID if --remote and all Organism objects created during
    the genomic context stage. Methods for generating summary tables and matrices for
    plotting are also stored in this class.

    Session objects can be dumped to/loaded from JSON to facilitate re-filtering/plotting.

    ::
        s = Session()
        with open("session.json", "w") as fp:
            s.to_json(fp)

        with open("session.json") as fp:
            s2 = Session.from_json(fp)

        s == s2
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

    def count_queries(self, hits, identity=False):
        """Count number of hits for each query in a collection of Hits.
        If identity=True, will return the maximum identity hit for each query.
        """

        def count(query):
            if identity:
                return max([h.identity if query == h.query else 0 for h in hits])
            return sum(query == h.query for h in hits)

        return [count(query) for query in self.queries]

    def form_matrices(self):
        """Form 2D count matrices for plotting."""

        def form_row(organism, accession, cluster):
            genus, species = organism.name.split(" ", 1)
            name = f"$\it{{{genus[0]}. {species}}}$ {organism.strain}"
            scaf = f"{accession}:{cluster[0].start}-{cluster[-1].end}"
            cnts = self.count_queries(cluster)
            idts = self.count_queries(cluster, identity=True)
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

        Parameters
        ----------
        organisms: list
        output: open file handle
        """
        return "\n\n\n".join(
            organism.summary(headers=headers, human=human)
            for organism in self.organisms
            if organism.count_hit_clusters() > 0
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
                [str(x) for x in self.count_queries(cluster, identity=identity)],
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
        """Tables, etc"""
        if form == "summary":
            table = self._summary(human=human, headers=headers, **kwargs)
        elif form == "binary":
            table = self._binary(human=human, headers=headers, **kwargs)
        else:
            raise ValueError("Expected 'summary' or 'binary'")
        print(table, file=fp)


class Organism(Serializer):
    """An organism.

    This class is used to represent a unique organism found during a cblaster search.
    Every strain (or lack thereof) is counted as being unique, and will be reported
    separately in cblaster results.
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
    """A genomic scaffold.

    This class represents a genomic scaffold belonging to an Organism object found
    during a cblaster search.
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
    """A BLAST hit.

    Stores hit scores and genomic context. It is first instantiated when parsing BLAST
    results, and is then updated with genomic coordinates after either querying the
    NCBI's IPG resource or a local JSON database.
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
        """Return a copy of this Hit object with any additional args."""
        copy = Hit(**self.__dict__)
        for key, val in kwargs.items():
            setattr(copy, key, val)
        return copy

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
