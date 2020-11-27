#!/usr/bin/env python3

""" this will act as an extension to the g2j package and allow parsing of embl files. Idealy
there the patterns in the original package are extended so only patterns have to be
"""

import re

import tempfile

import g2j
from g2j import genbank


PATTERNS = {
    "scaffold": re.compile(r"ID\s+?(?P<accession>\b[\w.-]+?);\s.+?\n//", re.DOTALL),
    "organism": re.compile(r"OS\s+?(?P<organism>\w.+?)[\n\r]"),
    "strain": re.compile(r'/strain=\"(?P<strain>[\w .-]+?)\"'),
    "features": re.compile(
        r"Key {13}Location/Qualifiers(.+?)(?:XX.+SQ {3}Sequence)", re.DOTALL
    ),
    "sequence": re.compile(r"SQ {3}Sequence[\w; ]+\s{6}(.*)", re.DOTALL),
    "identifier": re.compile(r'(protein_id|locus_tag|gene|ID)=\"([\w.:-]+?)\"'),
    "qualifier": re.compile(r'/(\w+?)=(?:\"(.+?)\"|(\d+?))', re.DOTALL),
}


def parse(query_file, feature_types=None, save_scaffold_sequence=True):
    """Parse an EMBL file

    Parameters:
        query_file (str): Path EMBL file
        feature_types (List): parsed features must be one of these types
        save_scaffold_sequence (Boolean): if the scaffold sequence should
        be saved
    Returns:
        organism (g2j.classes.Organims)
    """
    # has to be here or unexpected crashes occur
    genbank.PATTERNS = PATTERNS
    with open(query_file) as handle:
        new_text = handle.read().replace("FT   ", "     ")
    with tempfile.TemporaryFile(mode="w+") as tmp:
        tmp.writelines(new_text)
        tmp.seek(0)
        organism = genbank.parse(tmp, feature_types=feature_types, save_scaffold_sequence=save_scaffold_sequence)
    return organism
