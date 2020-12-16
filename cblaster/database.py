"""
This module handles creation of local JSON databases for non-NCBI lookups.
"""

import logging
import subprocess
import sqlite3
import warnings

from pathlib import Path
from multiprocessing import Pool

import gffutils
from gffutils import biopython_integration

from Bio import SeqIO, BiopythonParserWarning
from Bio.SeqFeature import FeatureLocation

from cblaster import helpers


# ignore malformed locus warnings
warnings.simplefilter('ignore', BiopythonParserWarning)


LOG = logging.getLogger("cblaster")


SCHEMA = """\
CREATE TABLE gene (
    id              INTEGER PRIMARY KEY,
    name 	    TEXT,
    start_pos       INTEGER,
    end_pos	    INTEGER,
    strand	    INTEGER,
    translation	    TEXT,
    scaffold        TEXT,
    organism        TEXT
);\
"""

QUERY = """\
SELECT
    id,
    name,
    start_pos,
    end_pos,
    strand,
    scaffold,
    organism
FROM
    gene
WHERE
    id IN ({})\
"""

FASTA = 'SELECT ">"||gene.id||"\n"||gene.translation||"\n" FROM gene'

INSERT = """\
INSERT INTO gene (
    name,
    start_pos,
    end_pos,
    strand,
    translation,
    scaffold,
    organism
)
VALUES
    (?, ?, ?, ?, ?, ?, ?)\
"""

GET_ORGANISM_IDS = "SELECT id, name FROM organism"

FASTA_SUFFIXES = (".fa", ".fsa", ".fna", ".fasta", ".faa")
GBK_SUFFIXES = (".gbk", ".gb", ".genbank", ".gbf", ".gbff")
GFF_SUFFIXES = (".gtf", ".gff", ".gff3")
EMBL_SUFFIXES = (".embl",)


def find_fasta(gff_path):
    """Finds a FASTA file corresponding to the given GFF path."""
    for suffix in FASTA_SUFFIXES:
        path = Path(gff_path).with_suffix(suffix)
        if path.exists():
            return path


def parse_fasta(path):
    with open(path) as fp:
        records = list(SeqIO.parse(fp, "fasta"))
    return records


def find_regions(directives):
    """Looks for ##sequence-region directives in a list of GFF3 directives."""
    regions = {}
    for directive in directives:
        if directive.startswith("sequence-region"):
            _, accession, start, end = directive.split(" ")
            regions[accession] = (int(start), int(end))
    return regions


def parse_gff(path):
    """Parses GFF and corresponding FASTA using GFFutils.

    - Parse gene/CDS features, translations & real coordinates
    - Convert to BioPython SeqRecords/SeqFeatures
    - All other formats just str -> SeqIO.parse(fp, 'filetype')
      then check for 'gff'
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

    # Find features for each record in the FASTA file
    for record in fasta:
        try:
            record_start, _ = regions[record.id]
        except KeyError:
            record_start = 1

        # Normalise Feature location based on ##sequence-region directive.
        # Necessary for extracted GFF3 files that still store coordinates
        # relative to the entire region. If no sequence-region directive
        # is found, assumes 1 (i.e. default sequence start).
        cds_features = []
        for feature in gff.region(seqid=record.id, featuretype=["gene", "CDS"]):
            feature = biopython_integration.to_seqfeature(feature)
            feature.location = FeatureLocation(
                feature.location.start - record_start - 1,
                feature.location.end - record_start,
                strand=feature.location.strand
            )
            if feature.type == "CDS":
                cds_features.append(feature)
            else:
                record.features.append(feature)

        if not record.features:
            raise ValueError(f"Found no CDS features in {record.id} [{path}]")

        # Merge CDS features into singular SeqFeature objects, add them to record
        previous = None
        for feature in sorted(cds_features, key=lambda f: f.location.start):
            seqid = feature.qualifiers["ID"][0]
            same_feature = previous == seqid
            if not previous:
                previous = seqid
            if same_feature:
                record.features[-1].location += feature.location
            else:
                record.features.append(feature)
                previous = seqid

        # Sort, then generate insertion tuples like with other formats
        record.features.sort(key=lambda f: f.location.start)

    return fasta


def parse_file(path):
    """Dispatches a given file path to the correct parser given its extension.

    Args:
        path (str): Path to genome file
    Returns:
        SeqRecord
    """
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix in GBK_SUFFIXES:
        file_type = "genbank"
    elif suffix in EMBL_SUFFIXES:
        file_type = "embl"
    elif suffix in GFF_SUFFIXES:
        records = parse_gff(path)
    else:
        raise ValueError(f"File {path} has invalid extension ({suffix})")
    with open(path) as fp:
        records = list(SeqIO.parse(fp, file_type))
    return dict(name=path.with_suffix("").name, records=records)


def init_sqlite_db(path, force=False):
    """Initialises a cblaster SQLite database file at a given path.

    Args:
        path (str): Path to write SQLite3 database
        force (bool): Overwrite pre-existing files at `path`
    """
    if Path(path).exists():
        if force:
            LOG.info("Overwriting pre-existing file at %s", path)
            Path(path).unlink()
        else:
            raise FileExistsError(f"File {path} already exists but force=False")
    else:
        LOG.info("Initialising cblaster SQLite3 database to %s", path)
    with sqlite3.connect(path) as con:
        con.executescript(SCHEMA)


def find_overlapping_location(feature, locations):
    """Finds the index of a gene location containing `feature`."""
    for index, (start, end) in enumerate(locations):
        if feature.location.start >= start and feature.location.end <= end:
            return index


def find_gene_name(qualifiers):
    """Finds a gene name in a dictionary of feature qualifiers."""
    for tag in ["locus_tag", "protein_id", "id", "gene", "name", "label"]:
        if tag in qualifiers:
            return qualifiers[tag]
    return "N.A."


def seqrecord_to_genes(record, source):
    """Generates insertion tuples for genes in a SeqRecord object.
    """
    features = [f for f in record.features if f.type == "CDS"]
    locations = [(f.location.start, f.location.end) for f in record.features if f.type == "gene"]
    genes = []
    for feature in features:
        qualifiers = {
            k: v[0] if isinstance(v, list) else v
            for k, v in feature.qualifiers.items()
        }
        name = find_gene_name(qualifiers)
        if "pseudo" in qualifiers:
            LOG.warning("%s is pseudogene, skipping", name)
            continue
        match = find_overlapping_location(feature, locations)
        if match:
            start, end = locations.pop(match)
        else:
            start, end = feature.location.start, feature.location.end

        # TODO: fix broken translations for GFF features
        translation = (
            qualifiers.pop("translation", None)
            or feature.extract(record.seq).translate()
        )

        if not translation:
            LOG.warning("Failed to find translation for %s, skipping", name)
            continue
        gene = (
            str(name),
            int(start),
            int(end),
            int(feature.location.strand),
            str(translation),
            str(record.id),  # scaffold accession
            str(source)  # organism name from source file name
        )
        genes.append(gene)
    return genes


def organisms_to_tuples(organisms):
    """Generates insertion tuples from a dictionary of parsed organisms."""
    tuples = []
    for organism in organisms:
        for record in organism["records"]:
            genes = seqrecord_to_genes(record, organism["name"])
            tuples.extend(genes)
    return tuples


def seqrecords_to_sqlite(tuples, database):
    """Writes a collection of SeqRecord objects to a cblaster SQLite database.
    """
    try:
        with sqlite3.connect(database) as con:
            cur = con.cursor()
            cur.executemany(INSERT, tuples)
    except sqlite3.IntegrityError:
        LOG.exception("Failed to insert %i records", len(tuples))


def sqlite_to_fasta(path, database):
    """Writes all proteins in `database` to `path` in FASTA format."""
    with sqlite3.connect(database) as con, open(path, "w") as fasta:
        cur = con.cursor()
        for (record,) in cur.execute(FASTA):
            fasta.write(record)


def query_database(ids, database):
    """Queries the cblaster SQLite3 database for a collection of gene IDs."""
    marks = ", ".join("?" for _ in ids)
    query = QUERY.format(marks)
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        return cur.execute(query, ids).fetchall()


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


def makedb(paths, database, force=False, cpus=None, batch=None):
    """Makedb module entry point"""
    LOG.info("Starting makedb module")

    if not (batch is None or isinstance(batch, int)):
        raise TypeError("batch should be None or int")
    if not (cpus is None or isinstance(cpus, int)):
        raise TypeError("cpus should be None or int")

    sqlite_path = Path(f"{database}.sqlite3")
    fasta_path = Path(f"{database}.fasta")
    dmnd_path = Path(f"{database}.dmnd")

    if sqlite_path.exists() or dmnd_path.exists():
        if force:
            LOG.info("Pre-existing files found, overwriting")
        else:
            raise RuntimeError("Existing files found but force=False")

    LOG.info("Initialising SQLite3 database at %s", sqlite_path)
    init_sqlite_db(sqlite_path, force=force)

    total_paths = len(paths)
    if batch is None:
        batch = total_paths
    path_groups = [paths[i: i + batch] for i in range(0, total_paths, batch)]
    LOG.info(
        "Parsing %i genome files, in %i batches of %i",
        total_paths,
        len(path_groups),
        batch
    )

    with Pool(cpus) as pool:
        for index, group in enumerate(path_groups, 1):
            LOG.info("Batch %i: %s", index, group)
            organisms = pool.map(parse_file, group)
            tuples = organisms_to_tuples(organisms)

            LOG.info("Saving %i genes", len(tuples))
            seqrecords_to_sqlite(tuples, sqlite_path)

    LOG.info("Writing FASTA to %s", fasta_path)
    sqlite_to_fasta(fasta_path, sqlite_path)

    LOG.info("Building DIAMOND database at %s", dmnd_path)
    diamond_makedb(fasta_path, dmnd_path)

    LOG.info("Done!")
