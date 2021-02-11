"""
This module handles creation of local JSON databases for non-NCBI lookups.
"""

import logging
import subprocess
import sqlite3

from pathlib import Path
from multiprocessing import Pool

from cblaster import helpers, sql
from cblaster import genome_parsers as gp


LOG = logging.getLogger("cblaster")


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
    with sqlite3.connect(path) as con:
        con.executescript(sql.SCHEMA)


def seqrecords_to_sqlite(tuples, database):
    """Writes a collection of SeqRecord objects to a cblaster SQLite database.

    Args:
        tuples (list): Gene insertion tuples
        database (str): Path to SQLite3 database
    """
    try:
        with sqlite3.connect(database) as con:
            cur = con.cursor()
            cur.executemany(sql.INSERT, tuples)
    except sqlite3.IntegrityError:
        LOG.exception("Failed to insert %i records", len(tuples))


def sqlite_to_fasta(path, database):
    """Writes all proteins in `database` to `path` in FASTA format.

    Args:
        path (str): Path to output FASTA file
        database (str): Path to SQLite3 database
    """
    with sqlite3.connect(database) as con, open(path, "w") as fasta:
        cur = con.cursor()
        for (record,) in cur.execute(sql.FASTA):
            fasta.write(record)


def _query(query, values, database, fetch="all"):
    with sqlite3.connect(database) as con:
        cur = con.cursor()
        query = cur.execute(query, values)
        return query.fetchall() if fetch == "all" else query.fetchone()


def query_sequences(ids, database):
    marks = ", ".join("?" for _ in ids)
    query = sql.SEQUENCE_QUERY.format(marks)
    return _query(query, ids, database)


def query_genes(ids, database):
    """Queries the cblaster SQLite3 database for a collection of gene IDs.

    Args:
        ids (list): Row IDs of genes being queried
        database (str): Path to SQLite3 database
    Returns:
        list: Result tuples returned by the query
    """
    marks = ", ".join("?" for _ in ids)
    query = sql.GENE_QUERY.format(marks)
    return _query(query, ids, database)


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
    return _query(query, [*names, scaffold, organism, start, end], database)


def query_nucleotides(scaffold, organism, start, end, database):
    """Queries a database for a """
    query = sql.SCAFFOLD_QUERY.format(start, end - start)
    return _query(query, [scaffold, organism], database, fetch="one")


def diamond_makedb(fasta, name):
    """Builds a DIAMOND database

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


def makedb(paths, database, force=False, cpus=None, batch=None):
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
    fasta_path = Path(f"{database}.fasta")
    dmnd_path = Path(f"{database}.dmnd")

    if sqlite_path.exists() or dmnd_path.exists():
        if force:
            LOG.info("Pre-existing files found, overwriting")
        else:
            raise RuntimeError("Existing files found but force=False")

    init_sqlite_db(sqlite_path, force=force)

    paths = gp.find_files(paths)
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

    with Pool(cpus) as pool:
        for index, group in enumerate(path_groups, 1):
            LOG.info("Batch %i: %s", index, [str(p) for p in group])
            organisms = pool.map(gp.parse_file, group)
            tuples = gp.organisms_to_tuples(organisms)

            LOG.info("Saving %i genes", len(tuples))
            seqrecords_to_sqlite(tuples, sqlite_path)

    LOG.info("Writing FASTA to %s", fasta_path)
    sqlite_to_fasta(fasta_path, sqlite_path)

    LOG.info("Building DIAMOND database at %s", dmnd_path)
    diamond_makedb(fasta_path, dmnd_path)

    LOG.info("Done!")
