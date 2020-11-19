#!/usr/bin/env python3

""" this will act as an extension to the g2j package and allow parsing of embl files. Idealy
there the patterns in the original package are extended so only patterns have to be
"""

import re


import g2j
from g2j import genbank
from g2j.classes import Organism, Scaffold, Feature


PATTERNS = {
    "scaffold": re.compile(r"ID\s+?(?P<accession>\b[\w.-]+?);\s.+?\n//", re.DOTALL),
    "organism": re.compile(r"OS\s+?(?P<organism>\w.+?)[\n\r]"),
    "strain": re.compile(r'/strain=\"(?P<strain>[\w .-]+?)\"'),
    "features": re.compile(
        r"Key {13}Location/Qualifiers(.+?)(?:SQ {3}Sequence)", re.DOTALL
    ),
    "sequence": re.compile(r"SQ {3}Sequence[\w; ]+\s{6}(.*)", re.DOTALL),
    "identifier": re.compile(r'(protein_id|locus_tag|gene|ID)=\"([\w.:-]+?)\"'),
    "qualifier": re.compile(r'/(\w+?)=(?:\"(.+?)\"|(\d+?))', re.DOTALL),
}

genbank.PATTERNS = PATTERNS


def parse_qualifier_block(text):
    """Parse qualifiers from qualifier block.

    Qualifiers are split by newline -> FT 19 spaces -> slash
    If value, Remove leading/trailing ", leading newline -> FT 19 spaces
    Otherwise it's boolean, e.g. /pseudo, so set True
    Store qualifiers of same type in lists, otherwise just str

    Parameters:
        text (string): string to extract the qualifiers from
    Returns:
        Dictionary of qualifiers with names of the qualifiers as keys
    """
    qualifiers = {}
    # TODO: with a pattern in the genbank file of genome2jason this function does not have to be overwritten
    gap = "\n" + "FT" + " " * 19
    for qualifier in text.split(f"{gap}/"):
        key, *value = qualifier.split("=")

        # Determine if boolean or str value
        if value:
            value = value[0].lstrip('"').strip('"\n').replace(gap, "")
        else:
            value = True

        # Save multiple qualifiers of same type in lists
        if key in qualifiers:
            if isinstance(qualifiers[key], list):
                qualifiers[key].append(value)
            else:
                qualifiers[key] = [qualifiers.pop(key), value]
        else:
            qualifiers[key] = value

    return qualifiers


def parse_feature_block(text, types=None):
    """Parse features from feature block.

    Block should resemble:
        ^     feature         1..100
        ^                     /qualifier1="value1"
        ^                     /qualifier2="value2"
        ^     feature2        200..500
        ...
    as matched using PATTERNS["features"].
    Separate feature blocks are identified by looking for differences
    in whitespace; new feature lines start with FT and 3 spaces, qualifiers
    with FT and 19.

    Parameters:
        text (string): the string the features have to be extracted from
        types (list): parsed features must be one of these types
    Returns:

    """
    # RegEx pattern for finding feature starts
    # TODO: with a pattern in the genbank file of genome2jason this function does not have to be overwritten
    pattern = re.compile(r"^FT {3}([\w0-9']+?)\s", re.M)
    features = []
    text_length = len(text)

    # Get positions of each unique sequence feature
    for match in pattern.finditer(text):
        feature = {
            "type": match.group(1),
            "interval": [match.end(), text_length],
        }
        if features:
            features[-1]["interval"][1] = match.start()
        features.append(feature)

    # Convert above to Feature objects
    for ix, feature in enumerate(features):

        # Get qualifier text corresponding to the current feature
        # The start+1 gets rid of the / in the first qualifier
        start, end = feature["interval"]
        block = text[start:end].split("/", 1)

        # Check for single line features
        if len(block) == 1:
            # remove the FT
            location = block[0].replace("FT", "")
            qualifiers = {}
        else:
            # remove the FT
            location = block[0].replace("FT", "")
            qualifiers = parse_qualifier_block(block[1])

        try:
            codon_start = int(qualifiers["codon_start"])
        except KeyError:
            codon_start = 1
        features[ix] = Feature(
            type=feature["type"],
            location=genbank.parse_location(location, codon_start),
            qualifiers=qualifiers,
        )

    if types:
        return [f for f in features if f.type in types]

    return features


def parse(handle, feature_types=None, save_scaffold_sequence=True):
    organism = Organism()
    for scaffold in genbank.scaffold_iter(handle.read()):

        text = scaffold.group(0)
        accession = scaffold.group("accession")

        if not organism.name:
            organism.name = genbank.find_pattern("organism", text)

        if not organism.strain:
            organism.strain = genbank.find_pattern("strain", text)

        feature_block = genbank.get_feature_block(text)

        if feature_block:
            scaffold = Scaffold(
                accession,
                genbank.get_scaffold_sequence(text) if save_scaffold_sequence else "",
                parse_feature_block(feature_block, types=feature_types),
            )
            organism.scaffolds.append(scaffold)
    return organism
