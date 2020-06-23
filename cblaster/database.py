"""
This module handles creation of local JSON databases for non-NCBI lookups.
"""

import json
import logging
import subprocess

from pathlib import Path

import g2j
from g2j import classes, genbank

from cblaster import helpers

LOG = logging.getLogger("cblaster")


class Database:
    """A cblaster database.

    This class acts as a wrapper around the genome2json classes (organism, scaffold,
    protein), consolidating them for use as a database.
    """

    def __init__(self, organisms=None):
        self.organisms = organisms if organisms else []

    def __iter__(self):
        return iter(self.organisms)

    def write_fasta(self, handle):
        """Formats organisms in the database to indexed FASTA format.
        Builds FASTA of each organism, then writes to given handle.
        """
        for i, organism in enumerate(self.organisms):
            fasta = ""
            for j, scaffold in enumerate(organism.scaffolds):
                for k, feature in enumerate(scaffold.features):
                    try:
                        sequence = feature.qualifiers["translation"]
                    except KeyError:
                        continue
                    fasta += f">{i}_{j}_{k}\n{sequence}\n"
            handle.write(fasta)

    @classmethod
    def from_files(cls, files):
        """Builds a new Database from a collection of GenBank files.

        For example:

        >>> db = Database.from_files(['path/to/file.gbk', 'path/to/file.gbk'])

        Only CDS features are parsed.
        """
        organisms = []
        LOG.info("Parsing %i files...", len(files))
        for index, file in enumerate(files, 1):
            with open(file) as handle:
                LOG.info("%i. %s", index, file)
                organism = g2j.genbank.parse(
                    handle,
                    feature_types=["CDS"],
                    save_scaffold_sequence=False
                )
                organisms.append(organism)
        return cls(organisms)

    def to_list(self):
        """Serialises all organisms in this database to dict."""
        return [organism.to_dict() for organism in self]

    def to_json(self, handle, indent=None):
        """Writes database to file in JSON format."""
        json.dump(self.to_list(), handle, indent=indent)

    @classmethod
    def from_json(cls, json_file):
        """Load a Database from JSON."""
        with open(json_file) as handle:
            js = json.load(handle)
        organisms = [g2j.classes.Organism.from_dict(d) for d in js]
        return cls(organisms)

    def makedb(self, name):
        """Convenience function to write FASTA and generate diamond DB"""
        fasta = f"{name}.faa"
        with open(fasta, "w") as handle:
            self.write_fasta(handle)
        diamond_makedb(fasta, name)


def diamond_makedb(fasta, name):
    """Builds a DIAMOND database from JSON.

    Args:
        fasta (str): Path to FASTA file containing protein sequences.
        name (str): Name for DIAMOND database.
    """
    diamond = helpers.get_program_path(["diamond", "diamond-aligner"])
    subprocess.run(
        [diamond, "makedb", "--in", fasta, "--db", name],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


def get_genbank_paths(folder):
    """Generates a collection of paths to GenBank files in a specified folder."""
    if not Path(folder).is_dir():
        raise ValueError("Expected valid folder")
    valid_extensions = (".gb", ".gbk", ".genbank")
    return [
        file
        for file
        in Path(folder).iterdir() if str(file).endswith(valid_extensions)
    ]
