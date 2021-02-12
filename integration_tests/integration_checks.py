"""Simple testing functions simply making sure that commands can run with certain values and some basic output checks

Warning: File comparissons and summary output rely heavily on diamond. These tests are made with diamond v 2.0.6
earlier or later versions will probably produce slightly different results and the provided databases will only
work with this exact version of diamond. That is why these versions are included"""

import os
import subprocess
from pathlib import Path
import platform
import pytest


CURRENT_DIR = None
TEST_FILE_DIR = None
COMPARISSON_FILE_DIR = None
OS_NAME = platform.system().lower()

POSSIBLE_FLAGS = ["silent"]


@pytest.fixture(scope="session", autouse=True)
def setup():
    global OS_NAME
    if OS_NAME == "darwin":
        raise SystemExit("Test for Mac OS are not yet configured.")

    global TEST_FILE_DIR, COMPARISSON_FILE_DIR, CURRENT_DIR
    # make sure the paths point the right way
    CURRENT_DIR = str(Path(__file__).resolve().parent)
    TEST_FILE_DIR = CURRENT_DIR + os.sep + "test_files"
    COMPARISSON_FILE_DIR = CURRENT_DIR + os.sep + "comparisson_files"
    os.chdir(CURRENT_DIR)

    # add the local versions of diamond to path
    old_path = os.environ["PATH"]
    os.environ["PATH"] = f"{str(Path('./diamond_files').resolve().absolute())}{os.pathsep}{old_path}"


def run_command(command):
    popen = subprocess.Popen(command, shell=True)
    stdout, stderr = popen.communicate()
    return popen.returncode, stderr, stdout


def confirm_files_present(files, directory):
    for file in os.listdir(directory):
        assert file in files


def compare_file_pairs(file_pairs, out_dir):
    for pair in file_pairs:
        compare_files(pair[0], pair[1], out_dir)


def compare_files(actual_file_path, expected_file_path, out_dir):
    with open(out_dir.join(actual_file_path), "r") as actual_file, \
            open(COMPARISSON_FILE_DIR + os.sep + expected_file_path, "r") as expected_file:
        for line_index, actual_expected in enumerate(zip(actual_file, expected_file)):
            actual_line, expected_line = actual_expected
            assert actual_line == expected_line


def test_gbk_query_local_search(tmpdir):
    command = \
        f"cblaster -d --testing search -m local -qf {TEST_FILE_DIR}{os.sep}test_query.gb -o " \
        f"{str(tmpdir.join('summary.txt'))} -db {TEST_FILE_DIR}{os.sep}test_database_{OS_NAME}.dmnd -ohh" \
        f" -ode , -odc 2 -osc -b {str(tmpdir.join('binary.txt'))} -bhh -bde _ -bdc 2 -bkey sum -bat coverage " \
        f" --blast_file {str(tmpdir.join('blast.txt'))} --ipg_file {str(tmpdir.join('ipgs.txt'))} " \
        f"-g 25000 -u 2 -mh 3 -r AEK75493.1 -me 0.01 -mi 30 -mc 50 -s {str(tmpdir.join('session.json'))} -ig -md 6000" \
        f" -mic 2 -p {str(tmpdir.join('plot.html'))}"
    return_code, sdterr, stdout = run_command(command)
    assert return_code == 0
    expected_files = {"binary.txt", "blast.txt", "ipgs.txt", "session.json", "plot.html", "summary.txt"}
    confirm_files_present(expected_files, tmpdir)
    compare_file_pairs([["summary.txt", "summary_local_gbk.txt"], ["binary.txt", "binary_local_gbk.txt"]], tmpdir)


def test_test():
    assert 0

#
# def search_local_commands():
#     global OUT_DIR, TEST_FILE_DIR
#
#     commands = [
#         # test gbk query in local mode with all options enabled
#         CommandTest(
#             f"..{os.sep}cblaster{os.sep}main.py -d search -m local -qf {TEST_FILE_DIR}{os.sep}test_query.gb -o "
#             f"{OUT_DIR}{os.sep}summary.txt -db {TEST_FILE_DIR}{os.sep}test_database_{OS_NAME}.dmnd -ohh"
#             f" -ode , -odc 2 -osc -b {OUT_DIR}{os.sep}binary.txt -bhh -bde _ -bdc 2 -bkey sum -bat coverage "
#             f" --blast_file {OUT_DIR}{os.sep}blast.txt --ipg_file {OUT_DIR}{os.sep}ipgs.txt "
#             f"-g 25000 -u 2 -mh 3 -r AEK75493.1 -me 0.01 -mi 30 -mc 50 -s {OUT_DIR}{os.sep}session.json -ig -md 6000"
#             f" -mic 2",
#             "test gbk query local",
#             [["summary.txt", "summary_local_gbk.txt"], ["binary.txt", "binary_local_gbk.txt"]]
#         ),
#         # test embl query in local mode with all options enabled
#         CommandTest(
#             f"cblaster -d search -m local -qf {TEST_FILE_DIR}{os.sep}test_query.embl -o "
#             f"{OUT_DIR}{os.sep}summary.txt -db {TEST_FILE_DIR}{os.sep}test_database_{OS_NAME}.dmnd -ohh"
#             f" -ode , -odc 2 -osc -b {OUT_DIR}{os.sep}binary.txt -bhh -bde _ -bdc 2 -bkey sum -bat coverage "
#             f" --blast_file {OUT_DIR}{os.sep}blast.txt --ipg_file {OUT_DIR}{os.sep}ipgs.txt "
#             f"-g 25000 -u 2 -mh 3 -r AEK75493.1 -me 0.01 -mi 30 -mc 50 -s {OUT_DIR}{os.sep}session.json -ig -md 6000"
#             f" -mic 25",
#             "test embl query local"
#         ),
#         # test fasta query in local mode with all options enabled
#         CommandTest(
#             f"cblaster -d search -m local -qf {TEST_FILE_DIR}{os.sep}test_query.fa -o "
#             f"{OUT_DIR}{os.sep}summary.txt -db {TEST_FILE_DIR}{os.sep}test_database_{OS_NAME}.dmnd -ohh"
#             f" -ode , -odc 2 -osc -b {OUT_DIR}{os.sep}binary.txt -bhh -bde _ -bdc 2 -bkey sum -bat coverage "
#             f" --blast_file {OUT_DIR}{os.sep}blast.txt --ipg_file {OUT_DIR}{os.sep}ipgs.txt "
#             f"-g 25000 -u 2 -mh 3 -r AEK75493.1 -me 0.01 -mi 30 -mc 50 -s {OUT_DIR}{os.sep}session.json",
#             "test fasta query local"
#         ),
#         # test query identifiers in local mode
#         CommandTest(
#             f"cblaster -d search -m local -qi AEK75490.1 AEK75490.1 AEK75500.1 AEK75516.1 AEK75516.1"
#             f" AEK75502.1 -o {OUT_DIR}{os.sep}summary.txt -db "
#             f"{TEST_FILE_DIR}{os.sep}test_database_{OS_NAME}.dmnd -ohh -ode , -odc 2 -osc -b"
#             f" {OUT_DIR}{os.sep}binary.txt -bhh -bde _ -bdc 2 -bkey sum -bat coverage "
#             f" --blast_file {OUT_DIR}{os.sep}blast.txt --ipg_file {OUT_DIR}{os.sep}ipgs.txt "
#             f"-g 25000 -u 2 -mh 3 -me 0.01 -mi 30 -mc 50 -s {OUT_DIR}{os.sep}session.json -ig -md 6000"
#             f" -mic 2",
#             "test query identifiers local"
#         ),
#         # test local session with all options enabled
#         CommandTest(
#             f"cblaster -d search -m local -qf {TEST_FILE_DIR}{os.sep}test_query.gb "
#             f"-s {TEST_FILE_DIR}{os.sep}test_session_local_embl_{OS_NAME}.json "
#             f"{TEST_FILE_DIR}{os.sep}test_session_local_gbk_{OS_NAME}.json -db"
#             f" {TEST_FILE_DIR}{os.sep}test_database_{OS_NAME}.dmnd -o {OUT_DIR}{os.sep}summary.txt -ig -md 6000"
#             f" -mic 2",
#             "test local session",
#             [["summary.txt", "summary_local_gbk_embl_combined.txt"]]
#         ),
#         # test recompute a local session
#         CommandTest(
#             f"cblaster -d search -m local -qf {TEST_FILE_DIR}{os.sep}test_query.gb "
#             f"-s {OUT_DIR}{os.sep}test_session_local_gbk_copy.json --recompute"
#             f" -db {TEST_FILE_DIR}{os.sep}test_database_{OS_NAME}.dmnd -o {OUT_DIR}{os.sep}summary.txt "
#             f"-g 50000 -u 5 -mh 3 -me 0.01 -mi 30 -mc 50 -b {OUT_DIR}{os.sep}binary.txt -ig -md 2500"
#             f" -mic 2",
#             "test recompute local",
#             [["summary.txt", "summary_local_gbk_recompute.txt"], ["binary.txt", "binary_local_gbk_recompute.txt"]]
#         )]
#     return commands
#
#
# def search_remote_commands():
#     global OUT_DIR, TEST_FILE_DIR
#
#     # load a remote session and make summary and binary files
#     commands = [
#         CommandTest(
#             f"cblaster -d search -m remote -qf {TEST_FILE_DIR}{os.sep}test_query.fa -s "
#             f"{TEST_FILE_DIR}{os.sep}test_session_remote_fa_{OS_NAME}.json"
#             f" -o {OUT_DIR}{os.sep}summary.txt -b {OUT_DIR}{os.sep}binary.txt -ig",
#             "load remote session",
#             [["summary.txt", "summary_remote_fa.txt"], ["binary.txt", "binary_remote_fa.txt"]]
#         ),
#         # recompute a remote session
#         CommandTest(
#             f"cblaster -d search -m remote -qf {TEST_FILE_DIR}{os.sep}test_query.fa -s "
#             f"{TEST_FILE_DIR}{os.sep}test_session_remote_fa_{OS_NAME}.json"
#             f" -o {OUT_DIR}{os.sep}summary.txt -b {OUT_DIR}{os.sep}binary.txt --recompute -g 50000 -u 7"
#             f" -mh 3 -me 0.01 -mi 20 -mc 60 --sort_clusters -ig -md 6000 -mic 3 -osc",
#             "recompute remote session",
#             [["summary.txt", "summary_remote_fa_recompute.txt"], ["binary.txt", "binary_remote_fa_recompute.txt"]]
#         )
#     ]
#     return commands
#
#
# def search_hmm_commands():
#     global OUT_DIR, TEST_FILE_DIR, CURRENT_DIR
#     if OS_NAME == "windows":
#         print("Skipping hmm tests because the hmmer program is not supported on windows systems.")
#         return []
#     commands = [
#         CommandTest(
#             f"cblaster search -m hmm -qp PF00698 -pfam {CURRENT_DIR}{os.sep} -db "
#             f"{TEST_FILE_DIR}{os.sep}test_database_{OS_NAME}.fasta -ohh -o {OUT_DIR}{os.sep}summary.txt"
#             f" -ode , -odc 2 -osc -b {OUT_DIR}{os.sep}binary.txt -bhh -bde _ -bdc 2 -bkey sum -bat coverage "
#             f"-g 25000 -u 2 -mh 3 -me 0.01 -s {OUT_DIR}{os.sep}session.json",
#             "hmm search"
#         ),
#         CommandTest(
#             f"cblaster search -m hmm -qp PF00698 -pfam {CURRENT_DIR}{os.sep} -db "
#             f"{TEST_FILE_DIR}{os.sep}test_database_{OS_NAME}.fasta -o {OUT_DIR}{os.sep}summary.txt"
#             f" -b {OUT_DIR}{os.sep}binary.txt -s {TEST_FILE_DIR}{os.sep}test_session_hmm_fa.json",
#             "load hmm session"
#         )
#     ]
#     return commands
#
#
# def search_combi_local_command():
#     global OUT_DIR, TEST_FILE_DIR, CURRENT_DIR
#     if OS_NAME == "windows":
#         print("Skipping hmm tests because the hmmer program is not supported on windows systems.")
#         return []
#     commands = [
#         CommandTest(
#             f"cblaster search -m combi_local -qp PF00698 -pfam {CURRENT_DIR}{os.sep} -qf"
#             f" {TEST_FILE_DIR}{os.sep}test_query.gb -db {TEST_FILE_DIR}{os.sep}test_database_{OS_NAME}.fasta "
#             f" {TEST_FILE_DIR}{os.sep}test_database_{OS_NAME}.dmnd -ohh -o {OUT_DIR}{os.sep}summary.txt"
#             f" -ode , -odc 2 -osc -b {OUT_DIR}{os.sep}binary.txt -bhh -bde _ -bdc 2 -bkey sum -bat coverage "
#             f"-g 25000 -u 2 -mh 3 -me 0.01 -s {OUT_DIR}{os.sep}session.json",
#             "combi_local search"
#         ),
#         CommandTest(
#             f"cblaster search -m combi_local -qp PF00698 -pfam {CURRENT_DIR}{os.sep} -db "
#             f"{TEST_FILE_DIR}{os.sep}test_database_{OS_NAME}.fasta -o {OUT_DIR}{os.sep}summary.txt"
#             f" -b {OUT_DIR}{os.sep}binary.txt -s {TEST_FILE_DIR}{os.sep}test_session_combi_local_fa.json",
#             "load combi local hmm session"
#         )
#     ]
#     return commands
#
#
# def makedb_commands():
#     global OUT_DIR, TEST_FILE_DIR
#     commands = [
#         # using embl and gbk files
#         CommandTest(
#             f"cblaster -d makedb {TEST_FILE_DIR}{os.sep}test_query.gb"
#             f" {TEST_FILE_DIR}{os.sep}test_query.embl -n {OUT_DIR}{os.sep}database -b 3 -c 100",
#             "makedb gbk, embl files",
#             [["database.fasta", "makedb_database_gbk_embl.fasta"]]
#         ),
#         # using gff files
#         CommandTest(
#             f"cblaster -d makedb {TEST_FILE_DIR}{os.sep}test_gff_v_maris.fna"
#             f" {TEST_FILE_DIR}{os.sep}test_gff_v_maris.gff -n {OUT_DIR}{os.sep}database -b 3 -c 100",
#             "makedb gff files"
#         )
#     ]
#     return commands
#
#
# def gne_commands():
#     global OUT_DIR, TEST_FILE_DIR
#     commands = [
#         CommandTest(
#             f"cblaster -d gne {TEST_FILE_DIR}{os.sep}test_session_local_gbk_{OS_NAME}.json --max_gap 200000 --samples 25"
#             f' --scale log -o {OUT_DIR}{os.sep}summary.txt -hh -d "\t" -e 2',
#             "gne with local gbk session",
#             [["summary.txt", "gne_local_summary.txt"]]
#         ),
#         CommandTest(
#             f"cblaster -d gne {TEST_FILE_DIR}{os.sep}test_session_remote_fa_{OS_NAME}.json --max_gap 250000 --samples 10"
#             f' --scale linear -o {OUT_DIR}{os.sep}summary.txt -hh -d "\t" -e 3',
#             "gne with remote fa session",
#             [["summary.txt", "gne_remote_summary.txt"]]
#         )
#     ]
#     return commands
#
#
# def extract_commands():
#     global OUT_DIR, TEST_FILE_DIR
#     commands = [
#         CommandTest(
#             f"cblaster -d extract {TEST_FILE_DIR}{os.sep}test_session_local_gbk_{OS_NAME}.json "
#             f"-or GCA_0002041.* -q AEK75517.1 AEK75496.1 AEK75490.1 -o"
#             f" {OUT_DIR}{os.sep}summary.txt -sc CP002638.1:2562030-2619476 -de _ -es -no",
#             "local session extraction",
#             [["summary.txt", "extract_local_summary.txt"]]
#         ),
#         CommandTest(
#             f"cblaster -d extract {TEST_FILE_DIR}{os.sep}test_session_remote_fa_{OS_NAME}.json -o "
#             f"{OUT_DIR}{os.sep}summary.txt -de : -es -or Verrucosispora.*",
#             "remote session extraction",
#             [["summary.txt", "extract_remote_summary.txt"]]
#         )
#     ]
#     return commands
#
#
# def extract_clusters_commands():
#     global OUT_DIR, TEST_FILE_DIR
#     commands = [
#         CommandTest(
#             f"cblaster -d extract_clusters {TEST_FILE_DIR}{os.sep}test_session_local_gbk_{OS_NAME}.json -c 1-2 4 -o"
#             f" {OUT_DIR} -pf test_ -f genbank -mc 5",
#             "local session cluster extraction",
#             [["test_cluster1.gbk", "extract_clusters_cluster1_local_genbank.gbk"]]
#         ),
#         CommandTest(
#             f"cblaster -d extract_clusters {TEST_FILE_DIR}{os.sep}test_session_remote_fa_{OS_NAME}.json -o {OUT_DIR} -or"
#             f" Verrucosispora.* -sc KF826676.1:5000-30000 -st 20 -f bigscape",
#             "remote session cluster extraction",
#             [["cluster11.gbk", "extract_clusters_cluster11_remote_bigscape.gbk"]]
#         )
#     ]
#     return commands
#
#
# def plot_clusters_commands():
#     # this check will open 2 plots that are deleter right after. There is unfortunately no way to stop this.
#     global OUT_DIR, TEST_FILE_DIR
#     commands = [
#         CommandTest(
#             f"cblaster -d plot_clusters {TEST_FILE_DIR}{os.sep}test_session_local_gbk_{OS_NAME}.json -c 1-2 4 -o"
#             f" {OUT_DIR}{os.sep}plot.html -mc 5",
#             "local session cluster plot"
#         ),
#         CommandTest(
#             f"cblaster -d plot_clusters {TEST_FILE_DIR}{os.sep}test_session_remote_fa_{OS_NAME}.json -o"
#             f" {OUT_DIR}{os.sep}plot.html -mc 5 -or Verrucosispora.* -sc KF826676.1:5000-30000 -st 20",
#             "remote session cluster plot"
#         )
#     ]
#     return commands


if __name__ == '__main__':
    pytest.main(args=['-s', os.path.abspath(__file__), "--log-cli-level", "DEBUG"])
#     # run tests of all the requested commands by specifying the command name separated by spaces
#     cmd_arguments = sys.argv[1:]
#     print("Running setup")
#     setup()
#     print("\nRunning tests:\n")
#     test_commands(cmd_arguments)
