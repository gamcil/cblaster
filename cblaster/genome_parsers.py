import warnings
import logging

from pathlib import Path

import gffutils
from gffutils import biopython_integration

from Bio import SeqIO, BiopythonParserWarning
from Bio.SeqFeature import FeatureLocation


# ignore malformed locus warnings
warnings.simplefilter('ignore', BiopythonParserWarning)

LOG = logging.getLogger("cblaster")

FASTA_SUFFIXES = (".fa", ".fsa", ".fna", ".fasta", ".faa")
GBK_SUFFIXES = (".gbk", ".gb", ".genbank", ".gbf", ".gbff")
GFF_SUFFIXES = (".gtf", ".gff", ".gff3")
EMBL_SUFFIXES = (".embl", ".emb")


def find_overlapping_location(feature, locations):
    """Finds the index of a gene location containing `feature`.

    Args:
        feature (SeqFeature): Feature being matched to a location
        locations (list): Start and end coordinates of gene features
    Returns:
        int: Index of matching start/end, if any
        None: No match found
    """
    for index, (start, end) in enumerate(locations):
        if feature.location.start >= start and feature.location.end <= end:
            return index


def find_gene_name(qualifiers):
    """Finds a gene name in a dictionary of feature qualifiers."""
    if not isinstance(qualifiers, dict):
        raise TypeError("Expected qualifier dictionary")
    for tag in ["locus_tag", "protein_id", "id", "gene", "name", "label"]:
        if tag in qualifiers:
            return qualifiers[tag]
    return "N.A."


def find_fasta(gff_path):
    """Finds a FASTA file corresponding to the given GFF path."""
    for suffix in FASTA_SUFFIXES:
        path = Path(gff_path).with_suffix(suffix)
        if path.exists():
            return path


def parse_fasta(path):
    with open(path) as fp:
        return list(SeqIO.parse(fp, "fasta"))


def find_regions(directives):
    """Looks for ##sequence-region directives in a list of GFF3 directives."""
    regions = {}
    for directive in directives:
        if directive.startswith("sequence-region"):
            _, accession, start, end = directive.split(" ")
            regions[accession] = (int(start), int(end))
    return regions


def parse_cds_features(features, record_start):
    cds = []
    gene = []
    for feature in features:
        feature = biopython_integration.to_seqfeature(feature)
        feature.location = FeatureLocation(
            feature.location.start - record_start,
            feature.location.end - record_start,
            strand=feature.location.strand
        )
        if feature.type == "CDS":
            cds.append(feature)
        else:
            gene.append(feature)
    return cds, gene


def merge_cds_features(features):
    # Merge CDS features into singular SeqFeature objects, add them to record
    features.sort(key=lambda f: f.location.start)
    merged, features = features[:1], features[1:]
    for feature in features:
        if merged[-1].qualifiers["ID"][0] == feature.qualifiers["ID"][0]:
            if feature.location.strand == 1:
                merged[-1].location += feature.location
            else:
                # Reverse strand locations must be in biological order
                old, new = merged[-1].location, feature.location
                merged[-1].location = new + old
        else:
            merged.append(feature)
    return merged


def parse_gff(path):
    """Parses GFF and corresponding FASTA using GFFutils.

    Args:
        path (str):
            Path to GFF file. Should have a corresponding FASTA file of the same
            name with a valid FASTA suffix (.fa, .fasta, .fsa, .fna, .faa).
    Returns:
        list: SeqRecord objects corresponding to each scaffold in the file
    """
    fasta = find_fasta(path)
    if not fasta:
        raise FileNotFoundError(f"Could not find partner FASTA file for {path}")

    # Parse FASTA and create GFFUtils database
    fasta = parse_fasta(fasta)
    gff = gffutils.create_db(
        str(path),
        ":memory:",
        force=True,
        merge_strategy="create_unique",
        sort_attribute_values=True
    )
    regions = find_regions(gff.directives)

    for record in fasta:
        # Normalise Feature location based on ##sequence-region directive.
        # Necessary for extracted GFF3 files that still store coordinates
        # relative to the entire region, not to the extracted FASTA.
        # If no sequence-region directive is found, assumes 1 (i.e. sequence start).
        cds, gene = parse_cds_features(
            gff.region(seqid=record.id, featuretype=["gene", "CDS"]),
            regions[record.id][0] - 1 if record.id in regions else 0
        )
        if not cds:
            LOG.warning("Found no CDS features in %s [%s]", record.id, path)
        record.features = sorted(
            [*gene, *merge_cds_features(cds)],
            key=lambda f: f.location.start
        )

    return fasta


def find_files(paths, recurse=True, level=0):
    files = []
    for path in paths:
        _path = Path(path)
        if _path.is_dir():
            if level == 0 or recurse:
                _files = find_files(_path.iterdir(), recurse=recurse, level=level + 1)
                files.extend(_files)
        else:
            ext = _path.suffix.lower()
            valid = ext in GBK_SUFFIXES + GFF_SUFFIXES + EMBL_SUFFIXES
            if _path.exists() and valid:
                files.append(_path)
    return files


def parse_file(path, to_tuples=False):
    """Dispatches a given file path to the correct parser given its extension.

    Args:
        path (str): Path to genome file
    Returns:
        dict: File name and list of SeqRecord objects corresponding to scaffolds in file
    """
    path = Path(path)
    name = path.with_suffix("").name
    suffix = path.suffix.lower()
    if suffix in GBK_SUFFIXES:
        file_type = "genbank"
    elif suffix in EMBL_SUFFIXES:
        file_type = "embl"
    elif suffix in GFF_SUFFIXES:
        return dict(name=name, records=parse_gff(path))
    elif suffix in FASTA_SUFFIXES:
        return dict(name=name, records=parse_fasta(path))
    else:
        raise ValueError(f"File {path} has invalid extension ({suffix})")
    with open(path) as fp:
        records = [
            seqrecord_to_tuples(record, path.stem) if to_tuples else record
            for record in SeqIO.parse(fp, file_type)
        ]
    return dict(name=name, records=records)


def iter_overlapping_features(features):
    if not features:
        return

    sorted_features = sorted(features, key=lambda f: f.location.start)
    first = sorted_features.pop(0)
    group, border = [first], first.location.end

    for feature in sorted_features:
        if feature.location.end <= border:
            group.append(feature)
            border = max(border, feature.location.end)
        else:
            yield group
            group, border = [feature], feature.location.end
    yield group


def seqrecord_to_tuples(record, source):
    features = [f for f in record.features if f.type in ("CDS", "gene")]

    rows = []
    for group in iter_overlapping_features(features):
        # Pull out gene and CDS features
        gene = [f for f in group if f.type == "gene"]
        cds = [f for f in group if f.type == "CDS"]

        # Get coordinates; prefer gene instead of CDS
        base = gene[0] if gene else cds[0]
        start = int(base.location.start)
        end = int(base.location.end)
        strand = int(base.location.strand) if base.location.strand else 1

        # Get name and translation, prefer CDS instead of gene
        name = find_gene_name(cds[0].qualifiers if cds else gene[0].qualifiers)[0]
        translation = (
            cds[0].qualifiers.pop("translation", None)[0]
            or str(cds[0].extract(record.seq).translate())
        ) if cds else ""

        # Keep track of record ID and source
        record_id = str(record.id)
        source_id = str(source)

        # Create the actual tuple
        row = ("gene", name, start, end, strand, translation, record_id, source_id)
        rows.append(row)

    scaffold = (
        "scaffold",
        str(record.id),
        0,
        len(record.seq),
        0,
        str(record.seq),
        str(record.id),
        str(source),
    )

    return [scaffold, *rows]


def organisms_to_tuples(organisms):
    """Generates insertion tuples from parsed organisms.

    Args:
        organisms (list): Organism dictionaries parsed by `parse_file`
    Returns:
        list: SQLite3 database insertion tuples for all genes
    """
    tuples = []
    for organism in organisms:
        for record in organism["records"]:
            genes = seqrecord_to_tuples(record, organism["name"])
            tuples.extend(genes)
    return tuples
