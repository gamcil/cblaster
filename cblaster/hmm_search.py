#!/usr/bin/env python3

"""
Hmmfetch and hmmsearch implementation
"""

import gzip
import subprocess
import logging

from datetime import datetime
from ftplib import FTP
from pathlib import Path

from Bio import SearchIO
from cblaster.classes import Hit
from cblaster import helpers

LOG = logging.getLogger(__name__)


def check_pfam_db(path):
    """Check if Pfam-A db exists else download

    Args:
        path: String, path where to check
    """
    path = Path(path)

    if path.exists() and not path.is_dir():
        raise FileExistsError("Expected directory")

    if not path.exists():
        path.mkdir()

    hmm = path / "Pfam-A.hmm.gz"
    dat = path / "Pfam-A.hmm.dat.gz"

    if hmm.exists() and dat.exists():
        LOG.info("Pfam database found")
    else:
        LOG.info("Downloading latest release from Pfam FTP")
        with FTP("ftp.ebi.ac.uk") as ftp:
            ftp.login()
            ftp.cwd("pub/databases/Pfam/current_release")
            ftp.retrbinary(f"RETR {hmm.name}", hmm.open("wb").write)
            ftp.retrbinary(f"RETR {dat.name}", dat.open("wb").write)

    return hmm, dat


def get_full_accession_number(dat_path, keys):
    """Get full accession number of Pfam profiles

    Args:
        keys: List, Strings of accession profiles numbers
        db_path: String, Path to dat.gz file with the full acc-nr
    Return:
        key_lines: List, string of full acc-number
    """
    valid_keys = []
    for line in gzip.open(dat_path, "rt"):
        if not line.startswith("#=GF AC"):
            continue
        *_, accession = line.strip().split(" ")
        for key in keys:
            if key in accession:
                valid_keys.append(accession)
    return valid_keys


def fetch_profiles(hmm, dat, keys):
    """Fetch hmm profiles from db and save in a file

    Args:
        db_path: String, path where db are stored
        keys_ls: String, Path to file with acc-nr
    Return:
        ls_keys: List, strings with acc-numbers
    """
    LOG.info("Fetching profiles from Pfam-A file")
    if not isinstance(hmm, Path):
        hmm = Path(hmm)

    LOG.info("Profiles found: %s", keys)
    profiles = []
    for key in keys:
        result = subprocess.run(["hmmfetch", hmm, key], stdout=subprocess.PIPE)
        profiles.append(result.stdout)

    output = hmm.parent / datetime.now().strftime("cblaster_%Y%m%d%H%M%S.hmm")
    with output.open("wb") as fp:
        for profile in profiles:
            fp.write(profile)

    LOG.info("Saved profiles to %s", output)
    return output


def run_hmmsearch(pfam, fasta, query):
    """Run the hmmsearch command

    Args:
        path_pfam: String, Path to the pfam database
        path_db: String, Path to db that will be searched for profiles
        ls_keys: List, string of pfam profile names
    Return:
        temp_res: List, String of result file names
    """
    LOG.info("Performing hmmsearch")
    output = query.with_suffix(".txt")
    try:
        subprocess.run(
            f"hmmsearch -o {output} {query} {fasta}",
            stdout=subprocess.PIPE,
            shell=True,
            check=True,
        )
    except subprocess.CalledProcessError:
        LOG.exception("hmmsearch failed!")
    return output


def parse_hmmer_output(results):
    """Parse hmmsearch output

    Args:
        file_list: List, string of file name of results that need parsing
    Return:
        hit_info: list of class objects, with information
                 - query, subject, identity, coverage, e-value, bit score
    """
    hit_info = []
    for record in SearchIO.parse(results, 'hmmer3-text'):
        if not record.hits:
            continue
        for hit in record.hits:
            hit_class = Hit(
                query=record.accession,  # Pfam id
                subject=hit.id,  # Hit id
                identity=None,  # Not present
                coverage=None,  # Not present
                evalue=hit.evalue,  # E-value of hit
                bitscore=hit.bitscore,  # Bit score of hit
            )
            hit_info.append(hit_class)
    if not hit_info:
        LOG.error("No hits have been found")
    return hit_info


def perform_hmmer(fasta, query_profiles, pfam, session):
    """Main of running a hmmer search

    Args:
        database_pfam: String, Path to pfam db
        database: String, Path to seqeunce db, in fasta or gbk format
        query_profiles: List, Pfam profiles needed to be searched
    Returns:
        hit_res: List of class objects with the hits

    """
    LOG.info("Starting hmmer search")

    # Make sure we can find hmmfetch and hmmsearch on PATH
    helpers.get_program_path(["hmmfetch", "hmmsearch"])

    # Find Pfam database (.dat and .hmm)
    hmm, dat = check_pfam_db(pfam)

    # Find real Pfam accessions
    session.queries = get_full_accession_number(dat, query_profiles)

    if not session.queries:
        LOG.error("No valid profiles could be selected")
        return

    # Extract HMM profiles from database
    query = fetch_profiles(hmm, dat, session.queries)

    # Run search
    results = run_hmmsearch(hmm, fasta, query)

    # Parse results and return
    return parse_hmmer_output(results)
