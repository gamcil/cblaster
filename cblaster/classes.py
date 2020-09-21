#!/usr/bin/env python3

"""
This module stores the classes (Organism, Scaffold, Hit) used in cblaster.
"""

import re
import json

from cblaster.formatters import (
    binary,
    summary,
    summarise_scaffold,
    summarise_organism,
)


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

    def __init__(self, queries=None, sequences=None, params=None, organisms=None):
        self.queries = queries if queries else []
        self.params = params if params else {}
        self.organisms = organisms if organisms else []
        self.sequences = sequences if sequences else {}

    def __add__(self, other):
        if not isinstance(other, Session):
            raise NotImplementedError("Expected Session object")
        if not self.queries == other.queries:
            raise ValueError("Query sequences do not match")
        return Session(
            queries=self.queries,
            sequences=self.sequences,
            params=self.params,
            organisms=self.organisms + other.organisms
        )

    def to_dict(self):
        return {
            "queries": self.queries,
            "sequences": self.sequences,
            "params": self.params,
            "organisms": [o.to_dict() for o in self.organisms],
        }

    @classmethod
    def from_file(cls, file):
        with open(file) as fp:
            s = cls.from_json(fp)
        return s

    @classmethod
    def from_files(cls, files):
        if len(files) == 1:
            return cls.from_file(files[0])
        first, *rest = files
        s = cls.from_file(first)
        for f in rest:
            s += cls.from_file(f)
        return s

    @classmethod
    def from_dict(cls, d):
        return cls(
            queries=d.get("queries", None),
            sequences=d.get("sequences", None),
            params=d.get("params", None),
            organisms=[Organism.from_dict(o) for o in d.get("organisms", [])],
        )

    def format(self, form, fp=None, **kwargs):
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
            table = summary(self, **kwargs)
        elif form == "binary":
            table = binary(self, **kwargs)
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
        total_subjects = sum(len(s.subjects) for s in self.scaffolds.values())
        return "ORGANISM: {} {} [{} subjects on {} scaffolds]".format(
            self.name, self.strain, total_subjects, total_scaffolds
        )

    @property
    def clusters(self):
        return [
            cluster
            for scaffold in self.scaffolds.values()
            for cluster in scaffold.clusters
        ]

    @property
    def total_hit_clusters(self):
        """Counts total amount of hit clusters in this Organism."""
        return sum(len(scaffold.clusters) for scaffold in self.scaffolds.values())

    def summary(self, decimals=4, hide_headers=True, delimiter=None):
        return summarise_organism(
            self,
            decimals=decimals,
            hide_headers=headers,
            delimiter=delimiter,
        )

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

    def __init__(self, accession, clusters=None, subjects=None):
        self.accession = accession
        self.subjects = subjects if subjects else []
        self.clusters = clusters if clusters else []

    def __str__(self):
        return "SCAFFOLD: {} [{} hits in {} clusters]".format(
            self.accession, len(self.subjects), len(self.clusters)
        )

    def summary(self, hide_headers=False, delimiter=None, decimals=4):
        return summarise_scaffold(
            self,
            decimals=decimals,
            hide_headers=headers,
            delimiter=delimiter,
        )

    def to_dict(self):
        return {
            "accession": self.accession,
            "subjects": [subject.to_dict() for subject in self.subjects],
            "clusters": [
                [self.subjects.index(subject) for subject in cluster]
                for cluster in self.clusters
            ],
        }

    @classmethod
    def from_dict(cls, d):
        subjects = [Subject.from_dict(subject) for subject in d["subjects"]]
        clusters = [[subjects[ix] for ix in cluster] for cluster in d["clusters"]]
        return cls(accession=d["accession"], subjects=subjects, clusters=clusters)


class Subject(Serializer):
    """A sequence representing one or more BLAST hits.

    This class is instantiated during the contextual lookup stage. It is
    important since it allows for subject sequences which hit >1 of
    the query sequences, while still staying non-redundant.

    Attributes:
        hits (list): Hit objects referencing this subject sequence.
        ipg (int): NCBI Identical Protein Group (IPG) id.
        start (int): Start of sequence on parent scaffold.
        end (int): End of sequence on parent scaffold.
        strand (str): Strandedness of the sequence ('+' or '-').
    """

    def __init__(self, hits=None, name=None, ipg=None, start=None, end=None, strand=None):
        self.hits = hits if hits else []
        self.ipg = ipg
        self.name = name
        self.start = int(start) if start is not None else None
        self.end = int(end) if end is not None else None
        self.strand = strand

    def __eq__(self, other):
        if not isinstance(other, Subject):
            raise NotImplementedError("Expected Subject object")
        return (
            set(self.hits) == set(other.hits)
            and self.ipg == other.ipg
            and self.start == other.start
            and self.end == other.end
            and self.strand == other.strand
        )

    def to_dict(self):
        return {
            "hits": [hit.to_dict() for hit in self.hits],
            "name": self.name,
            "ipg": self.ipg,
            "start": self.start,
            "end": self.end,
            "strand": self.strand
        }

    def values(self, decimals=4):
        records = []
        for hit in self.hits:
            record = (
                *hit.values(decimals),
                str(self.start),
                str(self.end),
                self.strand,
            )
            records.append(record)
        return records

    @classmethod
    def from_dict(cls, d):
        return cls(
            hits=[Hit.from_dict(h) for h in d["hits"]],
            name=d.get("name"),
            ipg=d.get("ipg"),
            start=d.get("start"),
            end=d.get("end"),
            strand=d.get("strand"),
        )


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
        bitscore
    ):
        self.query = query

        if "gb" in subject or "ref" in subject:
            subject = re.search(r"\|([A-Za-z0-9\._]+)\|", subject).group(1)

        self.subject = subject
        self.bitscore = float(bitscore)
        self.identity = float(identity)
        self.coverage = float(coverage)
        self.evalue = float(evalue)

    def __str__(self):
        return (
            f"Hit: {self.query} - {self.subject}:"
            f" {self.identity:.2%}/{self.coverage:.2%}"
        )

    def __key(self):
        return (
            self.query,
            self.bitscore,
            self.identity,
            self.coverage,
            self.evalue
        )

    def __hash__(self):
        return hash(self.__key())

    def __eq__(self, other):
        if not isinstance(other, Hit):
            raise NotImplementedError("Expected Hit object")
        return self.__key() == other.__key()

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
        ]

    def to_dict(self):
        return {
            "query": self.query,
            "subject": self.subject,
            "identity": self.identity,
            "coverage": self.coverage,
            "evalue": self.evalue,
            "bitscore": self.bitscore,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(**d)
