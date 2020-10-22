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

LOG = logging.getLogger(__name__)

def check_pfam_db(path):
    """Check f Pfam-A db exists else download

    :param path: String, path where to check
    :param file_names: list of strings, names of file in a list
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
                LOG.exception("Error: Path or file does not exists")
            except urllib.error.URLError or urllib.error.HTTPError:
                LOG.exception("Error: Internet connection problem")
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
    # Select the incomplete ones from the dat.gz info
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


def run_hmmsearch(profile_names, path, db_name=""):
    """Run the hmmsearch command

    :param profile_names: List, String of names of used profiles
    :param path: String, Path to folder where all db and profiles are
    :param db_name: String, Name of used database, needs to be in fasta.gz
                    format
    """
    LOG.info("Preforming hmmsearch")
    for prof in profile_names:
        command_run_hmmsearch = "hmmsearch -o {} {} {} ".format(prof +
                            "_results.txt", path + prof + ".hmm", path + db_name)
        os.subprocess.run(command_run_hmmsearch, shell=True)


def parse_hmmer_output(file=""):
    """Parse hmmsearch output

    :return: hit_info: Nested list, information about the hit results
                        - Hit_id, hit description, evalue, bit-score
    """
    hit_info = []
    for record in SearchIO.parse(file, 'hmmer3-text'):
        query_id = record.id
        hits = record.hits
        num_hits = len(hits)
        if num_hits > 0:
            for hit in hits:
                hit_id = hit.id  # hit sequence ID
                hit_description = hit.description  # hit sequence description
                current_evalue = hit.evalue  # hit-level e-value
                current_bitscore = hit.bitscore # hit-level score
                hit_info.append([query_id, hit_id, hit_description,
                                 current_evalue, current_bitscore])

    print(hit_info)
    return hit_info


def preform_hmmer(path_pfam=None,
                  path_db=None,
                  acc_profile=None):
    #1. Check if program exist else give error message and stop program
    if which("hmmfetch") is None or which("hmmsearch") is None:
        LOG.error("Hmmer could not be found in PATH")
        raise SystemExit

    LOG.info("Starting hmmer search")
    #2. run check_pfam_db
    check_pfam_db(path_pfam)

    #3. get_full_acc_number and run hmmfetch
    fetch_profiles(path_pfam, acc_profile)


    #4. run hmmsearch

    #5. Parse hmm output, needs to be the same as blast output
    # - query_id
    # - Subject id
    # - ide
