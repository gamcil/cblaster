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
from shutil import which
from cblaster.classes import Hit

LOG = logging.getLogger(__name__)

def check_pfam_db(path):
    """Check f Pfam-A db exists else download

    :param path: String, path where to check
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

    :param keys: List, Strings of accession profiles numbers
    :param db_path: String, Path to dat.gz file with the full acc-nr
    :return: key_lines: List, string of full acc-number
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

    :param db_path: String, path where db are stored
    :param keys_ls: String, Path to file with acc-nr
    :return ls_keys: List, strings with acc-numbers
    """
    LOG.info("Fetching profiles from Pfam-A file")
    ls_keys = get_full_accession_number(db_path, keys_ls)
    if not ls_keys:
        LOG.error("No valid profiles could be selected")
    else:
        for key in ls_keys:
            command_fetch_profile = "hmmfetch -o {} {} {}".format(db_path +
                              key + ".hmm", db_path + "Pfam-A.hmm.gz", key)
            subprocess.run(command_fetch_profile, stdout=subprocess.PIPE,
                           shell=True)
    LOG.info("Profiles found: %s", ls_keys)
    return ls_keys


def run_hmmsearch(path_pfam, path_db, ls_keys):
    """Run the hmmsearch command

    :param path_pfam: String, Path to the pfam database
    :param path_db: String, Path to db that will be searched for profiles
    :param ls_keys: List, string of pfam profile names
    :return: temp_res: List, String of result file names
    """
    LOG.info("Preforming hmmsearch")
    temp_res = []
    # TODO save results in a temp folder, delete when done
    for prof in ls_keys:
        command_run_hmmsearch = "hmmsearch -o {} {} {} ".format(prof +
                            "_results.txt", path_pfam + prof + ".hmm", path_db)
        subprocess.run(command_run_hmmsearch, stdout=subprocess.PIPE,
                       shell=True)
        temp_res.append(prof +"_results.txt")
    return temp_res


def parse_hmmer_output(file):
    """Parse hmmsearch output
    :param file: String, file name of results that need parsing
    :return: hit_info: Nested list, information about the hit results
                        - Hit_id, hit description, evalue, bit-score
    """
    hit_info = []
    for record in SearchIO.parse(file, 'hmmer3-text'):
        hits = record.hits
        num_hits = len(hits)
        if num_hits > 0:
            for hit in hits:
                ## TODO Parse correct ID from Results
                hit_class = Hit(
                    query=record.accession,  # Pfam id
                    subject=hit.id,  # Hit id
                    identity=None,  # Not present
                    coverage=None,  # Not present
                    evalue=hit.evalue,  # E-value of hit
                    bitscore=hit.bitscore,  # Bitscore of hit
                )
                hit_info.append(hit_class)
    if len(hit_info) == 0:
        LOG.error("No hits have been found")
    return hit_info


def preform_hmmer(path_pfam=None,
                  path_db=None,
                  acc_profile=None):
    #1. Check if program exist else give error message and stop program
    if which("hmmfetch") is None or which("hmmsearch") is None:
        LOG.error("Hmmer could not be found in PATH")

    LOG.info("Starting hmmer search")
    #2. run check_pfam_db
    check_pfam_db(path_pfam)

    #3. get_full_acc_number and run hmmfetch
    ls_keys = fetch_profiles(path_pfam, acc_profile)

    #4. run hmmsearch
    ls_res = run_hmmsearch(path_pfam, path_db, ls_keys)

    #5. Parse hmm output, needs to be the same as blast output
    ls_hits = []
    for res in ls_res:
        ls_hits.append(parse_hmmer_output(res))
    return ls_hits
