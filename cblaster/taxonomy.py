#!/usr/bin/env python3

"""
Generates taxonomic reports from NCBI taxonomy IDs
"""

import xml.etree.ElementTree as ET

from collections import Counter
from io import BytesIO

from Bio import Entrez

from cblaster.helpers import batch_function
from cblaster.classes import Serializer


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
    @batch_function(batch_size=100)
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
        return cls(d["taxa"], d["root"])
    
    def to_dict(self):
        return dict(taxa=self.taxa, root=self.root)