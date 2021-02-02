#!/usr/bin/env python3

"""
Hmmfetch and hmmsearch implementation
"""

import os
import gzip
import subprocess
import urllib.request
import urllib.error
import logging

from Bio import SearchIO
from cblaster.classes import Hit
from cblaster import helpers

LOG = logging.getLogger(__name__)


def check_pfam_db(path):
    """Check if Pfam-A db exists else download

    Args:
        path: String, path where to check
    """
    file_names = ["Pfam-A.hmm.gz", "Pfam-A.hmm.dat.gz"]
    url_ls = ["ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.1/Pfam-A.hmm.gz",
              "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.1/Pfam-A.hmm.dat.gz"]
    if os.path.exists(path + file_names[0]) and os.path.exists(path + file_names[1]):
        LOG.info("Pfam database found")
    else:
        LOG.info("Fetching database from Pfam release: 33.1 ")
        counter = 0
        for url in url_ls:
            try:
                urllib.request.urlretrieve(url, path + file_names[counter])
            except FileNotFoundError:
                LOG.error("Error: Path or file does not exists")
            except urllib.error.URLError or urllib.error.HTTPError:
                LOG.error("Error: Internet connection problem")
            counter += 1


def get_full_accession_number(db_path, keys):
    """Get full accession number of Pfam profiles

    Args:
        keys: List, Strings of accession profiles numbers
        db_path: String, Path to dat.gz file with the full acc-nr
    Return:
        key_lines: List, string of full acc-number
    """
    # Read dat.gz file with complete acc-numbers
    dat_gz_file = gzip.open(db_path + 'Pfam-A.hmm.dat.gz', 'r')
    content = str(dat_gz_file.read()).split("\\n")
    # Only appends to list when it is found in dat.gz file
    profile_ls = []
    for text in content:
        for key in keys:
            if key.strip() in text:
                profile_ls.append(text.split(" ")[-1])
    return profile_ls


def fetch_profiles(db_path, keys_ls):
    """Fetch hmm profiles from db and save in a file

    Args:
        db_path: String, path where db are stored
        keys_ls: String, Path to file with acc-nr
    Return:
        ls_keys: List, strings with acc-numbers
    """
    LOG.info("Fetching profiles from Pfam-A file")
    ls_keys = get_full_accession_number(db_path, keys_ls)
    if not ls_keys:
        LOG.error("No valid profiles could be selected")
        return "no matches"
    else:
        for key in ls_keys:
            command_fetch_profile = "hmmfetch -o {} {} {}".\
                format(db_path + key + ".hmm", db_path + "Pfam-A.hmm.gz", key)
            subprocess.run(command_fetch_profile, stdout=subprocess.PIPE,
                           shell=True)
    LOG.info("Profiles found: %s", ls_keys)
    return ls_keys


def run_hmmsearch(path_pfam, path_db, ls_keys):
    """Run the hmmsearch command

    Args:
        path_pfam: String, Path to the pfam database
        path_db: String, Path to db that will be searched for profiles
        ls_keys: List, string of pfam profile names
    Return:
        temp_res: List, String of result file names
    """
    LOG.info("Performing hmmsearch")
    temp_res = []
    for prof in ls_keys:
        command_run_hmmsearch = "hmmsearch -o {} {} {}".\
            format(prof + "_results.txt", path_pfam + prof + ".hmm", path_db)
        results = subprocess.run(command_run_hmmsearch, stdout=subprocess.PIPE, shell=True)
        if results.returncode != 0:
            LOG.error("hmmsearch command did not work")
        temp_res.append(prof + "_results.txt")
    return temp_res


def parse_hmmer_output(file_list):
    """Parse hmmsearch output

    Args:
        file_list: List, string of file name of results that need parsing
    Return:
        hit_info: list of class objects, with information
                 - query, subject, identity, coverage, e-value, bit score
    """
    hit_info = []
    for file in file_list:
        for record in SearchIO.parse(file, 'hmmer3-text'):
            hits = record.hits
            num_hits = len(hits)
            if num_hits > 0:
                for hit in hits:
                    hit_class = Hit(
                        query=record.accession,  # Pfam id
                        subject=hit.id,  # Hit id
                        identity=None,  # Not present
                        coverage=None,  # Not present
                        evalue=hit.evalue,  # E-value of hit
                        bitscore=hit.bitscore,  # Bit score of hit
                    )
                    hit_info.append(hit_class)
        if len(hit_info) == 0:
            LOG.error("No hits have been found")
    return hit_info


def perform_hmmer(
    database,
    query_profiles,
    database_pfam,
):
    """Main of running a hmmer search

    Args:
        database_pfam: String, Path to pfam db
        database: String, Path to seqeunce db, in fasta or gbk format
        query_profiles: List, Pfam profiles needed to be searched
    Returns:
        hit_res: List of class objects with the hits

    """
    # 1. Check if program exist else give error message and stop program
    helpers.get_program_path(["hmmfetch", "hmmsearch"])

    LOG.info("Starting hmmer search")
    # 2. run check_pfam_d
    check_pfam_db(database_pfam)

    # 3. get_full_acc_number and run hmmfetch
    ls_keys = fetch_profiles(database_pfam, query_profiles)
    if ls_keys == "no matches":
        return []

    # 4. run hmmsearch
    ls_res = run_hmmsearch(database_pfam, database, ls_keys)
    # 5. Parse hmm output, needs to be the same as blast output
    hit_res = parse_hmmer_output(ls_res)
    return hit_res
