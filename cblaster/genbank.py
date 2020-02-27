r"""
This module provides the parser for GenBank-format files.


The parser uses a series of regular expressions to extract information from a file.

Scaffold:
    LOCUS\s+?(?P<accession>\b[\w.-]+?)\s.+?//
        Find scaffold blocks in the file. This will match blocks of text, starting with
        LOCUS (first word of first line in a scaffold block), capturing anything until
        the scaffold terminates with //.

Organism:
    ORGANISM\s+?(?P<organism>\w[\w .-]+?)[\n\r]
        Find name of organism this DNA belongs to.

Strain:
    /strain="(?P<strain>([\w .-]+?))"
        Find strain name held in a /strain="STRAIN 123" type qualifier. Only check for
        this field specifically since there is no standard for storing in e.g. the
        DESCRIPTION or ORGANISM fields.

Protein:
    CDS\s+?
        Find start of CDS feature
    (.*?)
        Save any location modifiers (complement, join)
    (\d+?)\.\.
        Capture start position, e.g. 10..
    (?:.+?\.\.)*?
        Non-capture group, anything until .. before final
        number in location. Optional so that simple locations
        will also be parsed correctly.
    [<>]*?
        Account for potential truncation symbols
    (\d+?)
        Capture end position
    [<>)]*?
        Truncations, end parentheses...
    [\n\r]\s+?/
        Until next qualifier (starting with '/')
    (.*?)
        Capture all text until translation
    /translation="([A-Z\n\r\s ]+?)"
        Get translation

Identifier:
    (protein_id|locus_tag|gene|ID)="([\w.:-]+?)"
        Used with findall(), will retrieve any field matching these identifiers.
"""

import re

from pathlib import Path


PATTERNS = {
    "scaffold": re.compile(
        r"LOCUS\s+?(?P<accession>\b[\w.-]+?)\s.+?^//$", re.DOTALL|re.M
    ),
    "organism": re.compile(r"ORGANISM\s+?(?P<organism>\w[\w .-]+?)[\n\r]"),
    "strain": re.compile(r'/strain="(?P<strain>([\w .-]+?))"'),
    "protein": re.compile(
        r"CDS\s+?"
        r"([a-z<>(]*?)"  # join(complement( location modifiers, if any
        r"(\d+?)\.\."  # start
        r"(?:[0-9,.()]+?\.\.)*?"
        r"[<>)]*?"
        r"(\d+?)"  # end
        r"[<>)]*?"
        r"[\n\r]\s+?/"
        "(.*?)"  # capture everything, check for identifier in post
        r'/translation="([A-Z\n\r\s ]+?)"',  # translation
        re.DOTALL,
    ),
    "identifier": re.compile(r'(protein_id|locus_tag|gene|ID)="([\w.:-]+?)"'),
}


def scaffold_iter(text):
    """Thin wrapper around scaffold pattern."""
    yield from PATTERNS["scaffold"].finditer(text)


def find_pattern(pattern, text):
    """Find a match in a text block for a pattern in PATTERNS.

    If no match is found, as indicated by either IndexError or AttributeError being
    raised, this function will return None.
    """
    if pattern not in PATTERNS:
        raise ValueError("Invalid pattern specified")
    try:
        return PATTERNS[pattern].search(text).groups()[0]
    except (IndexError, AttributeError):
        return None


def find_best_identifier(text):
    """Find best protein identifier in a given text block.

    Coincidentally, reverse alphabetical order mirrors the preferred choice of
    protein identifier (i.e. protein_id -> locus_tag -> ID -> gene).

    So, reverse sort then take value of first match tuple.
    """

    ids = sorted(PATTERNS["identifier"].findall(text), reverse=True)

    if not ids:
        raise ValueError("No identifier could be found")

    return ids[0][1]


def find_strandedness(text):
    """Find strand of a feature given a location text block.

    Checks for 'complement', returning '-' if it is present, otherwise '+'.
    """
    return "-" if "complement" in text else "+"


def clean_sequence(text):
    """Clean translation sequence of whitespace and newline after found in regex."""
    return text.replace(" ", "").replace("\n", "")


def parse(handle):
    r"""Parse Proteins from a GenBank-format file.

    This function uses several regex patterns to extract gene positions. It first
    looks for scaffolds, by looking for blocks that start with LOCUS and end with //;
    each match is iterated.

    Another pattern is used inside each scaffold match to extract the start and end
    position of CDS features, as well as its protein ID and translation. Therefore, a
    minimal valid file should resemble something like:

    .. code-block::

        LOCUS     scaffold_1    ...
        ...
        ORGANISM  Organism name
        ...
        source    1..10000
                  /organism="Organism name"
                  /strain="STRAIN 123"
                  ...
        CDS       complement(join(1..100,200..300))
                  ...
                  /protein_id="Protein1"
                  ...
                  /translation="ABCDEF..."
        //
        ...

    If no organism or strain information can be identified, they will be set to 'NA',
    and CDS features will be parsed as normal.
    """

    organism, strain, scaffolds = None, None, []

    for scaffold in scaffold_iter(handle.read()):

        text = scaffold.group(0)
        accession = scaffold.group("accession")

        if not organism:
            organism = find_pattern("organism", text)

        if not strain:
            strain = find_pattern("strain", text)

        proteins = [
            {
                "index": index,
                "id": find_best_identifier(match[3]),
                "start": int(match[1]),
                "end": int(match[2]),
                "strand": find_strandedness(match[0]),
                "sequence": clean_sequence(match[-1]),
            }
            for index, match in enumerate(PATTERNS["protein"].findall(text))
        ]

        if proteins:
            scaffolds.append({"accession": accession, "proteins": proteins})

    if not organism:
        organism = "No_organism"

    if not strain:
        strain = "No_strain"

    return {
        "name": organism,
        "strain": strain,
        "file": handle.name,
        "scaffolds": scaffolds,
    }


def from_path(path):
    with open(path) as handle:
        organism = parse(handle)
    return organism


def get_genbank_paths(folder):
    """Generate a collection of paths to GenBank files in a specified folder."""
    if not Path(folder).is_dir():
        raise ValueError("Expected valid folder")
    valid_extensions = (".gb", ".gbk", ".genbank")
    return [
        file for file in Path(folder).iterdir() if str(file).endswith(valid_extensions)
    ]
