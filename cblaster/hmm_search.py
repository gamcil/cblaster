#!/usr/bin/env python3

"""
Hmmfetch and hmmsearch implementation
"""

import gzip
import subprocess
import logging
import re

from datetime import datetime
from ftplib import FTP
from pathlib import Path
from typing import Union, List, Collection, Set, Tuple

from Bio import SearchIO
from cblaster.classes import Hit, Session
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


def get_pfam_accession(
    dat_path: Union[str, Path],
    keys: Collection[str]
) -> Tuple[Set[str], Set[str]]:
    """Get full accession number of Pfam profiles

    Looks for keys in ID and AC attributes, such that accessions
    can be retrieved by name or accession.

    Args:
        keys: Strings of accession profiles numbers
        db_path: Path to dat.gz file with the full acc-nr
    Return:
        key_lines: List, string of full acc-number
    """
    keys = set(keys)
    valid_keys = set()
    name_attrs = ("#=GF ID", "#=GF AC")
    for line in gzip.open(dat_path, "rt"):
        if not line.startswith(name_attrs):
            continue
        *_, accession = line.strip().split(" ")
        if any(key in accession for key in keys if key not in valid_keys):
            valid_keys.add(accession)
    return valid_keys, keys.difference(valid_keys)


def read_profiles(files: Collection[str]) -> Collection[str]:
    """Reads in profile HMMs from a list of files."""
    profiles = [] 
    for file in files:
        with open(file) as fp:
            profile = fp.read()
            profiles.append(profile)
    return profiles


def get_profile_names(profiles: Collection[str]) -> Collection[str]:
    """Extracts names from profile HMMs using regular expressions."""
    names = []
    pattern = re.compile(r"^NAME\s+?(?P<name>\w.+?)$", re.M)
    for profile in profiles:
        matches = [match.group("name") for match in pattern.finditer(profile)]
        names.extend(matches)
    return names


def fetch_pfam_profiles(hmm, keys):
    """Fetch hmm profiles from db and save in a file

    Args:
        db_path: String, path where db are stored
        keys_ls: String, Path to file with acc-nr
    Return:
        ls_keys: List, strings with acc-numbers
    """
    if not isinstance(hmm, Path):
        hmm = Path(hmm)
    profiles = []
    for key in keys:
        result = subprocess.run(
            ["hmmfetch", hmm, key],
            stdout=subprocess.PIPE,
            encoding="utf-8",
        )
        profiles.append(result.stdout)
    return profiles


def write_profiles(profiles: Collection[str], output: str=None) -> str:
    """Writes a collection of profile HMMs to disk.

    If no output file is specified, will randomly generate a file name and save
    in the current working directory.

    Args:
        profiles: profile HMMs to write
        output: name of output file
    """
    if not output:
        output = datetime.now().strftime("cblaster_%Y%m%d%H%M%S.hmm")
    with open(output, "w") as fp:
        for profile in profiles:
            fp.write(profile)
    return output


def run_hmmsearch(fasta, query):
    """Run the hmmsearch command

    Args:
        path_pfam: String, Path to the pfam database
        path_db: String, Path to db that will be searched for profiles
        ls_keys: List, string of pfam profile names
    Return:
        temp_res: List, String of result file names
    """
    LOG.info("Performing hmmsearch")
    output = Path(query).with_suffix(".txt")
    try:
        subprocess.run(
            f"hmmsearch -o {output} {query} {fasta}",
            stdout=subprocess.PIPE,
            shell=True,
            check=True,
        )
    except subprocess.CalledProcessError:
        add_msg=""
        if fasta.endswith("gz"):
            add_msg = " Try running:  \n  gunzip {0} \nand then resubmit your cblaster command command, specifying:\n  -db {1}.\n\n".format(fasta, fasta.replace(".gz",""))
        LOG.exception("hmmsearch failed!"+add_msg)
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
                query=record.id,  # Pfam id
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


def group_profiles(profiles: Collection[str]) -> tuple:
    """Group input query profile HMMs by Pfam, custom or invalid.

    If the profile is found on disk, it will be loaded directly.
    If not found locally, but starts with PF, will try to extract from Pfam.
    Otherwise, marked as invalid.
    """
    local = []
    other = []
    for profile in profiles:
        if Path(profile).exists():
            local.append(profile)
        else:
            other.append(profile)
    return local, other


def perform_hmmer(
    fasta: str,
    query_profiles: List[str],
    pfam: str,
    session: Session,
) -> Union[Collection[Hit], None]:
    """Main of running a hmmer search

    Args:
        fasta: Path to database FASTA file
        query_profiles: Pfam names/accessions, or paths to profile HMM files
        pfam: Path to folder containing Pfam database
        session: cblaster search session
    Returns:
        List of class objects with the hits
    """
    LOG.info("Starting hmmer search")

    # Make sure we can find hmmfetch and hmmsearch on PATH
    helpers.get_program_path(["hmmfetch", "hmmsearch"])

    # Divide profiles into Pfam, local files and invalids
    local_profiles, other_profiles = group_profiles(query_profiles)

    # Read in/fetch all query profile HMMs.
    # Local profiles are just read in straight from file, Pfam accessions are
    # extracted from a local copy of the Pfam database (require .dat and .hmm).
    profiles = []
    if local_profiles:
        LOG.info("Loading local profiles: %s", local_profiles)
        profiles = read_profiles(local_profiles)
    if other_profiles:
        LOG.info("Fetching accessions from Pfam: %s", other_profiles)
        hmm, dat = check_pfam_db(pfam)
        accessions, invalid = get_pfam_accession(dat, other_profiles)
        if invalid:
            LOG.warning("Failed to fetch profiles from Pfam: %s", invalid)
        pfam_profiles = fetch_pfam_profiles(hmm, accessions)
        profiles.extend(pfam_profiles)
    if not profiles:
        LOG.error("No valid profiles could be selected")
        return
    query = write_profiles(profiles)
    LOG.info("Profiles written to: %s", query)

    # Save query profile HMM names
    session.queries = get_profile_names(profiles)

    # Run search
    results = run_hmmsearch(fasta, query)

    # Parse results and return
    return parse_hmmer_output(results)
