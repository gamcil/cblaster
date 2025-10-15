"""
This module handles creation of local JSON databases for non-NCBI lookups.
"""

import gzip
import logging
import subprocess
import sqlite3 as SQLITE
import functools

from pathlib import Path
from multiprocessing import Pool

from cblaster import helpers, sql
from cblaster import genome_parsers as gp

LOG = logging.getLogger(__name__)

try:
    import genomicsqlite as SQLITE
except (ImportError, ModuleNotFoundError, SQLITE.OperationalError) as exc:
    LOG.warning("Importing genomicsqlite failed, falling back to SQLite3")


def init_sqlite_db(path, force=False):
    """Initialises a cblaster SQLite3 database file at a given path.

    Args:
        path: Path to write SQLite3 database
        force: Overwrite pre-existing files at `path`

    Raises:
        FileExistsError: If `path` already exists but `force` is False
    """
    if Path(path).exists():
        if force:
            LOG.info("Overwriting pre-existing file at %s", path)
            Path(path).unlink()
        else:
            raise FileExistsError(f"File {path} already exists but force=False")
    else:
        LOG.info("Initialising cblaster SQLite3 database to %s", path)
    with SQLITE.connect(str(path)) as con:
        con.executescript(sql.SCHEMA)
        con.commit()


def seqrecords_to_sqlite(tuples, database):
    """Writes a collection of SeqRecord objects to a cblaster SQLite database.

    Args:
        tuples (list): Gene insertion tuples
        database (str): Path to SQLite3 database
    """
    try:
        with SQLITE.connect(str(database)) as con:
            cur = con.cursor()
            cur.executemany(sql.INSERT, tuples)
            con.commit()
        con.close()
    except SQLITE.IntegrityError:
        LOG.exception("Failed to insert %i records", len(tuples))


def sqlite_to_fasta(path, database, compress):
    """Writes all proteins in `database` to `path` in FASTA format.

    Args:
        path (str): Path to output FASTA file
        database (str): Path to SQLite3 database
    """
    if compress:
        handler = gzip.open(path, "wt")
    else:
        handler = open(path, 'w')

    with SQLITE.connect(str(database)) as con, handler  as fasta:
        cur = con.cursor()
        for (record,) in cur.execute(sql.FASTA):
            fasta.write(record)


def _query(query, database, values=None, fetch="all"):
    with SQLITE.connect(str(database)) as con:
        cur = con.cursor()
        query = cur.execute(query, values) if values else cur.execute(query)
        return query.fetchall() if fetch == "all" else query.fetchone()


def query_sequences(ids, database):
    inner = ", ".join(str(idx) for idx in ids)
    query = sql.SEQUENCE_QUERY.format(inner)
    return _query(query, database)


def query_genes(ids, database):
    """Queries the cblaster SQLite3 database for a collection of gene IDs.

    Args:
        ids (list): Row IDs of genes being queried
        database (str): Path to SQLite3 database
    Returns:
        list: Result tuples returned by the query
    """
    inner = ", ".join(str(idx) for idx in ids)
    query = sql.GENE_QUERY.format(inner)
    return _query(query, database)


def query_intermediate_genes(
    names, start, end, scaffold, organism, database, local=False
):
    """Queries the cblaster SQLite3 database for a collection of intermediate genes.

    These are the genes between start and stop that are not part of the names list

    Args:
        names (list): a list of names that are part of one cluster
        start (int): the minimal start a gene can have to be considered intermediate
        end (int): the maximum end a gene can have to be considered intermediate
        database (str): Path to SQLite3 database
    Returns:
        list: Result tuples returned by the query
    """
    marks = ", ".join("?" for _ in names)
    query = sql.INTERMEDIATE_GENES_QUERY.format(marks)
    return _query(query, database, values=[*names, scaffold, organism, start, end])


def query_nucleotides(scaffold, organism, start, end, database):
    """Queries a database for a """
    query = sql.SCAFFOLD_QUERY.format(start, end - start)
    return _query(query, database, values=[scaffold, organism], fetch="one")


def diamond_makedb(fasta, name, cpus):
    """Builds a DIAMOND database

    Args:
        fasta (str): Path to FASTA file containing protein sequences.
        name (str): Name for DIAMOND database.
    """
    diamond = helpers.get_program_path(["diamond", "diamond-aligner"])
    subprocess.run(
        [diamond, "makedb", "--in", str(fasta), "--db", name, "--threads", str(cpus)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


def makedb(paths, database, force=False, cpus=None, batch=None, compress=False):
    """makedb module entry point.

    Will parse genome files in `paths` and create:

        1. `database`.sqlite3
        SQLite3 database used for looking up genome context of hit genes

        2. `database`.dmnd
        DIAMOND search database

        3. `database`.fasta
        FASTA file containing all protein sequences in parsed genomes

    Args:
        paths (list): Paths to genome files to build database from
        database (str): Base name for database files
        force (bool): Overwrite pre-existing database files
        cpus (int):
            Number of CPUs to use when parsing genome files.
            By default, all available cores will be used.
        batch (int):
            Number of organisms to parse at once before saving to database.
            Helpful when dealing with larger/many genome files.
    """
    LOG.info("Starting makedb module")

    if not (batch is None or isinstance(batch, int)):
        raise TypeError("batch should be None or int")
    if not (cpus is None or isinstance(cpus, int)):
        raise TypeError("cpus should be None or int")

    sqlite_path = Path(f"{database}.sqlite3")
    dmnd_path = Path(f"{database}.dmnd")

    if compress:
        fasta_ext = '.fasta.gz'
    else:
        fasta_ext = '.fasta'

    fasta_path = Path(f"{database}{fasta_ext}")

    if sqlite_path.exists() or dmnd_path.exists():
        if force:
            LOG.info("Pre-existing files found, overwriting")
        else:
            raise RuntimeError("Existing files found but force=False")

    init_sqlite_db(sqlite_path, force=force)

    paths = gp.find_files(paths)
    if len(paths) == 0:
        raise RuntimeError("No valid files provided; expected genbank, embl or gff with accompanying fasta file (can be gzipped). Alternatively, provide one .txt file with one file per line.")

    # If a text file was provided, read and parse
    if len(paths)==1 and Path(paths[0]).suffix == ".txt":
        with open(paths[0]) as f:
            lines = f.read().splitlines()
        paths = gp.find_files(lines)

    total_paths = len(paths)
    if batch is None:
        batch = total_paths
    path_groups = [paths[i : i + batch] for i in range(0, total_paths, batch)]

    LOG.info(
        "Parsing %i genome files, in %i batches of %i",
        total_paths,
        len(path_groups),
        batch,
    )
    try:
        func = functools.partial(gp.parse_file, to_tuples=True)
        with Pool(cpus) as pool:
            for index, group in enumerate(path_groups, 1):
                LOG.info("Processing batch %i", index)
                for file in group:
                    LOG.info("  %s", file.name)
                tuples = []
                for organism in pool.imap(func, group):
                    for records in organism["records"]:
                        tuples.extend(records)
                LOG.info("Saving %i genes", len(tuples))
                seqrecords_to_sqlite(tuples, sqlite_path)
    except Exception:
        LOG.error("File parsing failed, exiting...", exc_info=True)

    LOG.info("Writing FASTA to %s", fasta_path)
    sqlite_to_fasta(fasta_path, sqlite_path, compress)

    LOG.info("Building DIAMOND database at %s", dmnd_path)
    diamond_makedb(fasta_path, dmnd_path, cpus)

    LOG.info("Done!")

