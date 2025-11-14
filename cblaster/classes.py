#!/usr/bin/env python3

"""
This module stores the classes (Organism, Scaffold, Hit) used in cblaster.
"""

import re
import json
import itertools
import xml.etree.ElementTree as ET

from abc import ABC, abstractmethod
from collections import Counter
from io import BytesIO

from Bio import Entrez

from cblaster.helpers import batch_function
from cblaster.formatters import (
    binary,
    summary,
    summarise_scaffold,
    summarise_organism,
)

class Serializer(ABC):
    """JSON serialisation mixin class.

    Classes that inherit from this class should implement `to_dict` and
    `from_dict` methods.
    """

    @abstractmethod
    def to_dict(self):
        """Serialises class to dict."""
        raise NotImplementedError

    @classmethod
    @abstractmethod
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
        sequences (dict): Query sequence translations
        query (Cluster): cblaster Cluster object for query
    """

    def __init__(
        self, queries=None, sequences=None, params=None, organisms=None, query=None, taxonomy=None,
    ):
        self.queries = queries if queries else []
        self.params = params if params else {}
        self.organisms = organisms if organisms else []
        self.sequences = sequences if sequences else {}
        self.taxonomy = taxonomy if taxonomy else {}
        self.query = query

    def __add__(self, other):
        if not isinstance(other, Session):
            raise NotImplementedError("Expected Session object")
        if not self.queries == other.queries:
            raise ValueError("Query sequences do not match")
        return Session(
            queries=self.queries,
            query=self.query,
            sequences=self.sequences,
            taxonomy=self.taxonomy,
            params=self.params,
            organisms=self.organisms + other.organisms,
        )
        
    def build_taxonomy(self):
        self.taxonomy = Taxonomy.from_session(self)

    def to_dict(self):
        return {
            "query": self.query.to_dict(save_subjects=True),
            "queries": self.queries,
            "params": self.params,
            "organisms": [o.to_dict() for o in self.organisms],
            "taxonomy": self.taxonomy.to_dict() if self.taxonomy else {}
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
            query=Cluster.from_dict(d["query"]),
            queries=d.get("queries", None),
            params=d.get("params", None),
            organisms=[Organism.from_dict(o) for o in d.get("organisms", [])],
            taxonomy=Taxonomy.from_dict(d.get("taxonomy", {}))
        )

    def format(self, form, fp=None, **kwargs):
        """Generates a summary table.

        Args:
            form (str): Type of table to generate ('summary' or 'binary').
            fp (file handle): File handle to write to.
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

    def __init__(self, name, strain, tax_id=None, scaffolds=None):
        self.name = name
        self.tax_id = tax_id
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
            self, decimals=decimals, hide_headers=hide_headers, delimiter=delimiter,
        )

    @property
    def full_name(self):
        """The full name (including strain) of the organism.
        Note: if strain found in name, returns just name.
        """
        if not self.name:
            name = "No organism"
            return f"{name} {self.strain}" if self.strain else name
        else:
            return f"{self.name} {self.strain}" if self.strain else self.name

    def to_dict(self):
        return {
            "name": self.name,
            "tax_id": self.tax_id,
            "strain": self.strain,
            "scaffolds": [scaffold.to_dict() for scaffold in self.scaffolds.values()],
        }

    @classmethod
    def from_dict(cls, d):
        return cls(
            d["name"],
            d["strain"],
            tax_id=d["tax_id"],
            scaffolds={
                scaffold["accession"]: Scaffold.from_dict(scaffold)
                for scaffold in d["scaffolds"]
            },
        )


class Scaffold(Serializer):
    """A genomic scaffold containing hits found in a cblaster search.

    Attributes:
        accession (str): Name of this scaffold, typically NCBI accession.
        subjects (list): Subject objects located on this scaffold.
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

    def add_clusters(self, subject_lists, query_sequence_order=None):
        """Add clusters to this scaffold

        After clusters are added they are sorted based on score

        Args:
            subject_lists (list): a list of lists of Subject objects that are
            form a clusters
            query_sequence_order (list): list of sequences of the order in the query file, is
            only provided if the query has a meningfull order (gbk, embl files).
        """
        for subjects in subject_lists:
            indices = [self.subjects.index(subject) for subject in subjects]
            cluster = Cluster(
                indices,
                subjects,
                query_sequence_order=query_sequence_order,
            )
            self.clusters.append(cluster)
        self.clusters.sort(key=lambda x: x.score, reverse=True)

    def remove_subject(self, subject):
        """Safely remove a subject from a cluster by removing it from the cluster as well.

        Args:
            subject (Subject): cblaster Subject object
        """
        remove_index = self.subjects.index(subject)
        for cluster in self.clusters:
            cluster.remove_subject(subject, remove_index)
            if len(cluster.subjects) == 0:
                self.clusters.remove(cluster)
        del self.subjects[remove_index]

    def summary(self, hide_headers=False, delimiter=None, decimals=4):
        return summarise_scaffold(
            self, decimals=decimals, hide_headers=hide_headers, delimiter=delimiter,
        )

    def to_dict(self):
        return {
            "accession": self.accession,
            "subjects": [subject.to_dict() for subject in self.subjects],
            "clusters": [cluster.to_dict() for cluster in self.clusters],
        }

    @classmethod
    def from_dict(cls, d):
        subjects = [Subject.from_dict(subject) for subject in d["subjects"]]
        clusters = [None] * len(d["clusters"])
        for index, cluster in enumerate(d["clusters"]):
            cluster_subjects = [subjects[ix] for ix in cluster["indices"]]
            clusters[index] = Cluster.from_dict(cluster, subjects=cluster_subjects)
        return cls(accession=d["accession"], subjects=subjects, clusters=clusters)


class Cluster(Serializer):
    """A cluster of subjects on the same scaffold

    Attributes:
        indices (list): indexes of the subjects in the list of subjects
        of the parent scaffold
        subjects (list): Subject objects that are in this cluster. Note:
        These are not serialised for this cluster
        intermediate_genes (list):
        start (int): The start coordinate of the cluster on the parent scaffold
        end (int): The end coordinate of the cluster on the parent scaffold
        number (int): number that is unique for each cluster in order to identify them
    """

    NUMBER = itertools.count()

    def __init__(
        self,
        indices=None,
        subjects=None,
        intermediate_genes=None,
        query_sequence_order=None,
        score=None,
        start=None,
        end=None,
        number=None,
    ):
        self.indices = indices if indices else []
        self.subjects = subjects if subjects else []
        self.intermediate_genes = intermediate_genes if intermediate_genes else []
        self.score = score if score else self.calculate_score(query_sequence_order)
        self.start = start
        self.end = end
        self.number = number if number is not None else next(self.NUMBER)

        if self.subjects and not (self.start or self.end):
            self.start = self.subjects[0].start
            self.end = self.subjects[-1].end

    def __iter__(self):
        return iter(self.subjects)

    def __len__(self):
        return len(self.subjects)

    def __eq__(self, other):
        if not isinstance(other, Cluster):
            raise NotImplementedError("Expected Cluster object")
        return (set(self.subjects) == set(other.subjects)
                and self.score == other.score)

    def __hash__(self):
        return hash(id(self))

    @property
    def sequences(self):
        return {
            subject.name: subject.sequence
            for subject in self.subjects
        }

    @property
    def names(self):
        return [subject.name for subject in self.subjects]

    def remove_subject(self, subject, scaffold_index):
        """Safely remove a subject from a cluster.

        This is important when subjects become empty when recomputing a session with different treshholds

        Args:
            subject (Subject): cblaster Subject object
            scaffold_index (int): the index of the subject in the scaffold it is saved in
        """
        try:
            remove_index = self.subjects.index(subject)
        except ValueError:
            return
        for index, index_value in enumerate(self.indices):
            if index_value > scaffold_index:
                self.indices[index] -= 1
        del self.indices[remove_index]
        del self.subjects[remove_index]

    def __calculate_synteny_score(self, query_sequence_order):
        if not query_sequence_order:
            return 0
        positions = []
        for index, subject in enumerate(self.subjects):
            best_hit = max(subject.hits, key=lambda h: h.bitscore)
            pair = (query_sequence_order.index(best_hit.query), index)
            positions.append(pair)
        score = 0
        for index, position in enumerate(positions[:-1]):
            query, subject = position
            next_query, next_subject = positions[index + 1]
            if abs(query - next_query) < 2 and abs(query - next_query) == abs(
                subject - next_subject
            ):
                score += 1
        return score

    def __calculate_bitscore(self):
        return sum(
            max(hit.bitscore for hit in subject.hits)
            for subject in self.subjects
            if subject.hits
        )

    @property
    def intermediate_start(self):
        """The start of the cluster taking the intermediate genes into account"""
        if not self.intermediate_genes:
            return self.start
        return min(*[s.start for s in self.intermediate_genes], self.start)

    @property
    def intermediate_end(self):
        """The end of the cluster taking the intermediate genes into account"""
        if not self.intermediate_genes:
            return self.end
        return max(*[s.end for s in self.intermediate_genes], self.end)

    def calculate_score(self, query_sequence_order=None):
        """Calculate the score of the current cluster

        The score is based on accumulated blastbitscore, total amount of hits against the
        query and a synteny score if query sequence order is provided. If there are multiple
        hits in a subject the hit with the top bitscore is selected for the calculation.

        Args:
            query_sequence_order (list): list of sequences of the order in the query file, is
            only provided if the query has a meningfull order (gbk, embl files).
        Returns:
            a float
        """
        synteny_score = self.__calculate_synteny_score(query_sequence_order)
        bitscore = self.__calculate_bitscore()
        return bitscore / 10000 + len(self.subjects) + synteny_score

    def to_dict(self, save_subjects=False):
        return {
            "indices": self.indices,
            "subjects": [s.to_dict() for s in self.subjects] if save_subjects else [],
            "intermediate_genes": [g.to_dict() for g in self.intermediate_genes],
            "score": self.score,
            "start": self.start,
            "end": self.end,
            "number": self.number,
        }

    @classmethod
    def from_dict(cls, d, subjects=None):
        return cls(
            indices=d["indices"],
            subjects=subjects or [Subject.from_dict(d) for d in d["subjects"]],
            intermediate_genes=[Subject.from_dict(d) for d in d["intermediate_genes"]],
            score=d["score"],
            start=d["start"],
            end=d["end"],
            number=d["number"],
        )


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

    def __init__(
        self,
        id=None,
        hits=None,
        name=None,
        ipg=None,
        start=None,
        end=None,
        strand=None,
        sequence=None,
    ):
        self.id = id
        self.hits = hits if hits else []
        self.ipg = ipg
        self.name = name
        self.start = int(start) if start is not None else None
        self.end = int(end) if end is not None else None
        self.strand = strand
        self.sequence = sequence

    def __key(self):
        # make equals behaviour of higher classes consistent with different instances
        return (
            self.id,
            tuple(sorted(self.hits, key=lambda x: x.bitscore)),
            self.ipg,
            self.start,
            self.end,
            self.strand
        )

    def __eq__(self, other):
        if not isinstance(other, Subject):
            raise NotImplementedError("Expected Subject object")
        return self.__key() == other.__key()

    def __hash__(self):
        return hash(self.__key())

    def to_dict(self):
        return {
            "id": self.id,
            "hits": [hit.to_dict() for hit in self.hits],
            "name": self.name,
            "ipg": self.ipg,
            "start": self.start,
            "end": self.end,
            "strand": self.strand,
            "sequence": self.sequence,
        }

    def values(self, decimals=4):
        start = str(self.start)
        end = str(self.end)
        strand = "+" if self.strand == 1 else "-"
        # If hits, this was found during the search 
        # Otherwise, it's an intermediate
        if len(self.hits) > 0:
            return [
                (*hit.values(decimals), start, end, strand)
                for hit in self.hits
            ]
        record = (
            "intermediate",
            self.name,
            "-", "-",
            "-", "-",
            start,
            end,
            strand,
        )
        return [record]

    @classmethod
    def from_dict(cls, d):
        return cls(
            id=d.get("id"),
            hits=[Hit.from_dict(h) for h in d["hits"]],
            name=d.get("name"),
            ipg=d.get("ipg"),
            start=d.get("start"),
            end=d.get("end"),
            strand=d.get("strand"),
            sequence=d.get("sequence"),
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
    """

    def __init__(self, query, subject, identity, coverage, evalue, bitscore):
        self.query = query

        if "gb" in subject or "ref" in subject:
            subject = re.search(r"\|([A-Za-z0-9\._]+)\|", subject).group(1)
        # Made id & Coverage a None type, hmmer does not have those values
        self.subject = subject
        self.bitscore = float(bitscore)
        self.identity = float(identity) if identity is not None else None
        self.coverage = float(coverage) if coverage is not None else None
        self.evalue = float(evalue)

    def __str__(self):
        return (
            f"Hit: {self.query} - {self.subject}:"
            f" {self.identity}/{self.coverage:}"
        )

    def __key(self):
        return self.query, self.bitscore, self.identity, self.coverage, self.evalue

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
            f"{round(self.identity, decimals):g}" if self.identity is not None else "",
            f"{round(self.coverage, decimals):g}" if self.coverage is not None else "",
            f"{self.evalue:.{decimals}g}" if self.evalue is not None else "",
            f"{round(self.bitscore, decimals):g}" if self.bitscore is not None else "",
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


class Taxonomy(Serializer):
    def __init__(self, taxa=None, root=None):
        self.taxa = taxa if taxa else {}
        self.root = root if root else None
        
    @staticmethod
    def parse_taxonomy_xml(records):
        """Parses Entrez.read() record"""
        entries = [
            dict(
                taxid=int(tax["TaxId"]),
                parent=int(tax["ParentTaxId"]),
                rank=tax.get("Rank", "no rank"),
                name=tax["ScientificName"],
                lineage=[
                    (int(x["TaxId"]), x["ScientificName"], x.get("Rank", "no rank"))
                    for x in tax.get("LineageEx", [])
                ],
                children=[]
            )
            for tax in records
        ]
        return entries

    @staticmethod
    def parse_xml_file(tax_xml_file):
        with open(tax_xml_file) as fp:
            raw_data = fp.read()
            xml_data = Entrez.read(BytesIO(raw_data.encode()))
            dicts = Taxonomy.parse_taxonomy_xml(xml_data)
        return dicts

    @staticmethod
    def fetch_taxonomy_dicts(batch):
        handle = Entrez.efetch(db="taxonomy", id=",".join(map(str, batch)), retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        entries = Taxonomy.parse_taxonomy_xml(records)
        return entries

    def build_taxonomy_tree(self, dicts):
        """Builds initial tree from list of taxonomy entry dictionaries"""
        self.taxa = { taxon["taxid"] : taxon for taxon in dicts}
        for tid, taxon in self.taxa.items():
            parent = taxon["parent"]
            if parent in self.taxa:
                self.taxa[parent]['children'].append(tid)
                
    def fill_missing_taxa(self):
        """Adds in intermediary nodes from lineages of fetched records"""
        for record in list(self.taxa.values()):
            lineage = record["lineage"]
            for i in range(len(lineage) - 1, -1, -1):
                lineage_id, lineage_name, lineage_rank = lineage[i]
                if lineage_id not in self.taxa:
                    self.taxa[lineage_id] = dict(
                        taxid=lineage_id,
                        parent=None,
                        rank=lineage_rank,
                        name=lineage_name,
                        lineage=[],
                        children=[],
                        count=0
                    )
                if i > 0:
                    parent_id, parent_name, parent_rank = lineage[i - 1]
                    if parent_id not in self.taxa: 
                        self.taxa[parent_id] = dict(
                            taxid=parent_id,
                            parent=None,
                            rank=parent_rank,
                            name=parent_name,
                            lineage=[],
                            children=[],
                            count=0
                        )
                    self.taxa[lineage_id]["parent"] = parent_id
                    if lineage_id not in self.taxa[parent_id]["children"]:
                        self.taxa[parent_id]["children"].append(lineage_id)
            if lineage:
                parent_id = lineage[-1][0]
                record["parent"] = parent_id
                if record["taxid"] not in self.taxa[parent_id]["children"]:
                    self.taxa[parent_id]["children"].append(record["taxid"])
                    
    def set_root(self):
        """Finds and saves the root node in the hierarchy"""
        self.root = next(iter(self.taxa))  # Pick any taxid (e.g., first in the dict)
        while "parent" in self.taxa[self.root] and self.taxa[self.root]["parent"] in self.taxa:
            self.root = self.taxa[self.root]["parent"]

    def set_counts(self, counts):
        """Sets initial count values from a Counter()"""
        for tid, count in counts.items():
            self.taxa[tid]["count"] = count
            
    def accumulate_counts(self, tid):
        """Computes cumulative counts for all nodes in the hierarchy"""
        def _accumulate(tid):
            node = self.taxa[tid]
            total = node.get("count", 0)
            for child_id in node.get("children", []):
                total += _accumulate(child_id)
            node["accumulated"] = total
            return total
        _accumulate(tid)

    def build(self, session, tax_xml_file=None):
        """Builds taxonomy from a Session"""
        counts = Counter(o.tax_id for o in session.organisms for _ in o.clusters)
        tids = list(counts)
        if (tax_xml_file):
            dicts = Taxonomy.parse_xml_file(tax_xml_file)
        else:
            dicts = Taxonomy.fetch_taxonomy_dicts(tids)
        self.build_taxonomy_tree(dicts)
        self.fill_missing_taxa()
        self.set_root()
        self.set_counts(counts)
        self.accumulate_counts(self.root)
    
    def format_kraken2(self):
        """Generates Kraken2-style taxonomy report"""

        root = self.taxa[self.root]
        root_total = root.get("accumulated", root.get("count", 0))

        def _format(tid, depth=0):
            node = self.taxa[tid]
            acc = node.get("accumulated", node.get("count", 0))
            raw = node.get("count", 0)
            percent_raw = 100 * raw / root_total if root_total else 0
            indent = "  " * depth
            line = f"{percent_raw:.2f}\t{acc}\t{raw}\t{node['rank'].lower()}\t{tid}\t{indent}{node['name']}"
            lines = [line]
            for child_id in sorted(node["children"], key=lambda x: self.taxa[x]["name"]):
                lines += _format(child_id, depth + 1)
            return lines
        
        report = [
            f"0.00\t0\t0\tno rank\t0\tunclassified",
            f"100.00\t{root_total}\t0\tno rank\t1\troot",
            *_format(self.root)
        ]
        return '\n'.join(report)
    
    @classmethod
    def from_session(cls, session, tax_xml_file=None):
        t = cls()
        t.build(session, tax_xml_file=tax_xml_file)
        return t
    
    @classmethod
    def from_dict(cls, d):
        return cls(d.get("taxa", {}), d.get("root", None))
    
    def to_dict(self):
        return dict(taxa=self.taxa, root=self.root)