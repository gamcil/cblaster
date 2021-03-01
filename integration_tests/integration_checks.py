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


def confirm_files_present(expected_files, directory):
    presen_files = os.listdir(directory)
    for file in expected_files:
        assert file in presen_files


def compare_file_pairs(file_pairs, out_dir):
    for pair in file_pairs:
        compare_files(pair[0], pair[1], out_dir)


def compare_files(actual_file_path, expected_file_path, out_dir):
    with open(out_dir.join(actual_file_path), "r") as actual_file, \
            open(COMPARISSON_FILE_DIR + os.sep + expected_file_path, "r") as expected_file:
        try:
            assert actual_file.read() == expected_file.read()
        # make the error a bit more readable by telling what files are not matching.
        except AssertionError as e:
            raise type(e)(str(e) + f" {actual_file.name} and {expected_file.name} dont match.")


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
    expected_files = {"binary.txt", "blast.txt", "session.json", "plot.html", "summary.txt"}
    confirm_files_present(expected_files, tmpdir)
    compare_file_pairs([["summary.txt", "summary_local_gbk.txt"], ["binary.txt", "binary_local_gbk.txt"]], tmpdir)


def test_embl_query_local_search(tmpdir):
    command = \
        f"cblaster -d --testing search -m local -qf {TEST_FILE_DIR}{os.sep}test_query.embl -o " \
        f"{str(tmpdir.join('summary.txt'))} -db {TEST_FILE_DIR}{os.sep}test_database_{OS_NAME}.dmnd -ohh" \
        f" -ode , -odc 2 -osc -b {str(tmpdir.join('binary.txt'))} -bhh -bde _ -bdc 2 -bkey sum -bat coverage " \
        f" --blast_file {str(tmpdir.join('blast.txt'))} --ipg_file {str(tmpdir.join('ipgs.txt'))} " \
        f"-g 25000 -u 2 -mh 3 -r AEK75493.1 -me 0.01 -mi 30 -mc 50 -s {str(tmpdir.join('session.json'))} -ig -md 6000" \
        f" -mic 2 -p {str(tmpdir.join('plot.html'))}"
    return_code, sdterr, stdout = run_command(command)
    assert return_code == 0
    expected_files = {"binary.txt", "blast.txt", "session.json", "plot.html", "summary.txt"}
    confirm_files_present(expected_files, tmpdir)


def test_fasta_query_local_search(tmpdir):
    command = \
        f"cblaster -d --testing search -m local -qf {TEST_FILE_DIR}{os.sep}test_query.fa -o " \
        f"{str(tmpdir.join('summary.txt'))} -db {TEST_FILE_DIR}{os.sep}test_database_{OS_NAME}.dmnd -ohh" \
        f" -ode , -odc 2 -osc -b {str(tmpdir.join('binary.txt'))} -bhh -bde _ -bdc 2 -bkey sum -bat coverage " \
        f" --blast_file {str(tmpdir.join('blast.txt'))} --ipg_file {str(tmpdir.join('ipgs.txt'))} " \
        f"-g 25000 -u 2 -mh 3 -r AEK75493.1 -me 0.01 -mi 30 -mc 50 -s {str(tmpdir.join('session.json'))} " \
        f"-p {str(tmpdir.join('plot.html'))}"
    return_code, sdterr, stdout = run_command(command)
    assert return_code == 0
    expected_files = {"binary.txt", "blast.txt", "session.json", "plot.html", "summary.txt"}
    confirm_files_present(expected_files, tmpdir)


def test_query_ids_local_search(tmpdir):
    command = \
        f"cblaster -d --testing search -m local -qi AEK75490.1 AEK75490.1 AEK75500.1 AEK75516.1 AEK75516.1 -o " \
        f"{str(tmpdir.join('summary.txt'))} -db {TEST_FILE_DIR}{os.sep}test_database_{OS_NAME}.dmnd -ohh" \
        f" -ode , -odc 2 -osc -b {str(tmpdir.join('binary.txt'))} -bhh -bde _ -bdc 2 -bkey sum -bat coverage " \
        f" --blast_file {str(tmpdir.join('blast.txt'))} --ipg_file {str(tmpdir.join('ipgs.txt'))} " \
        f"-g 25000 -u 2 -mh 3 -r AEK75493.1 -me 0.01 -mi 30 -mc 50 -s {str(tmpdir.join('session.json'))} -ig -md 6000" \
        f" -mic 2 -p {str(tmpdir.join('plot.html'))}"
    return_code, sdterr, stdout = run_command(command)
    assert return_code == 0
    expected_files = {"binary.txt", "blast.txt", "session.json", "plot.html", "summary.txt"}
    confirm_files_present(expected_files, tmpdir)


def test_load_session_local_search(tmpdir):
    command = \
        f"cblaster -d search -m local -qf {TEST_FILE_DIR}{os.sep}test_query.gb "\
        f"-s {TEST_FILE_DIR}{os.sep}test_session_local_embl_{OS_NAME}.json " \
        f"{TEST_FILE_DIR}{os.sep}test_session_local_gbk_{OS_NAME}.json -db" \
        f" {TEST_FILE_DIR}{os.sep}test_database_{OS_NAME}.dmnd -o {str(tmpdir.join('summary.txt'))} -ig -md 6000" \
        f" -mic 2"
    return_code, sdterr, stdout = run_command(command)
    assert return_code == 0
    expected_files = {"summary.txt"}
    confirm_files_present(expected_files, tmpdir)
    compare_file_pairs([["summary.txt", "summary_local_gbk_embl_combined.txt"]], tmpdir)


def test_recompute_session_local_search(tmpdir):
    command = \
        f"cblaster -d --testing search -m local -qf {TEST_FILE_DIR}{os.sep}test_query.gb " \
        f"-s {str(tmpdir.join('session.json'))} --recompute" \
        f" -db {TEST_FILE_DIR}{os.sep}test_database_{OS_NAME}.dmnd -o {str(tmpdir.join('summary.txt'))} " \
        f"-g 50000 -u 5 -mh 3 -me 0.01 -mi 30 -mc 50 -b {str(tmpdir.join('binary.txt'))} -ig -md 2500" \
        f" -mic 2 -p {str(tmpdir.join('plot.html'))}"
    return_code, sdterr, stdout = run_command(command)
    assert return_code == 0
    expected_files = {"binary.txt", "session.json", "plot.html", "summary.txt"}
    confirm_files_present(expected_files, tmpdir)
    compare_file_pairs([["summary.txt", "summary_local_gbk_recompute.txt"],
                        ["binary.txt", "binary_local_gbk_recompute.txt"]], tmpdir)


def test_load_remote_session(tmpdir):
    command = \
        f"cblaster -d --testing search -m remote -qf {TEST_FILE_DIR}{os.sep}test_query.fa -s " \
        f"{TEST_FILE_DIR}{os.sep}test_session_remote_fa_{OS_NAME}.json" \
        f" -o {str(tmpdir.join('summary.txt'))} -b {str(tmpdir.join('binary.txt'))} -ig" \
        f" -p {str(tmpdir.join('plot.html'))}"
    return_code, sdterr, stdout = run_command(command)
    assert return_code == 0
    expected_files = {"binary.txt", "plot.html", "summary.txt"}
    confirm_files_present(expected_files, tmpdir)
    compare_file_pairs([["summary.txt", "summary_remote_fa.txt"],
                        ["binary.txt", "binary_remote_fa.txt"]], tmpdir)


def test_recompute_remote_session(tmpdir):
    command = \
        f"cblaster -d --testing search -m remote -qf {TEST_FILE_DIR}{os.sep}test_query.fa -s " \
        f"{TEST_FILE_DIR}{os.sep}test_session_remote_fa_{OS_NAME}.json" \
        f" -o {str(tmpdir.join('summary.txt'))} -b {str(tmpdir.join('binary.txt'))} --recompute -g 50000 -u 7" \
        f" -mh 3 -me 0.01 -mi 20 -mc 60 --sort_clusters -ig -md 6000 -mic 3 -osc  -p {str(tmpdir.join('plot.html'))}"
    return_code, sdterr, stdout = run_command(command)
    assert return_code == 0
    expected_files = {"binary.txt", "plot.html", "summary.txt"}
    confirm_files_present(expected_files, tmpdir)
    compare_file_pairs([["summary.txt", "summary_remote_fa_recompute.txt"],
                        ["binary.txt", "binary_remote_fa_recompute.txt"]], tmpdir)


def test_search_hmm(tmpdir):
    command = \
        f"cblaster search -m hmm -qp PF00491 PF05593 -db {TEST_FILE_DIR}{os.sep}test_database_{OS_NAME}.fasta " \
        f"- pfam {TEST_FILE_DIR} -o {str(tmpdir.join('summary.txt'))} -s {str(tmpdir.join('session.json'))} " \
        f"-b {str(tmpdir.join('binary.txt'))} -p {str(tmpdir.join('plot.html'))}"
    return_code, sdterr, stdout = run_command(command)
    assert return_code == 0
    expected_files = {"binary.txt", "session.json", "plot.html", "summary.txt"}
    confirm_files_present(expected_files, tmpdir)


def test_embl_gbk_makedb(tmpdir):
    command = \
        f"cblaster -d --testing makedb {TEST_FILE_DIR}{os.sep}test_query.gb" \
        f" {TEST_FILE_DIR}{os.sep}test_query.embl -n {str(tmpdir.join('database'))} -b 3 -c 100"
    return_code, sdterr, stdout = run_command(command)
    assert return_code == 0
    expected_files = {"database.dmnd", "database.fasta", "database.sqlite3"}
    confirm_files_present(expected_files, tmpdir)
    compare_file_pairs([["database.fasta",  "makedb_database_gbk_embl.fasta"]], tmpdir)


def test_gff_fasta_makedb(tmpdir):
    command = \
        f"cblaster -d --testing makedb {TEST_FILE_DIR}{os.sep}test_gff_v_maris.fna" \
        f" {TEST_FILE_DIR}{os.sep}test_gff_v_maris.gff -n {str(tmpdir.join('database'))} -b 3 -c 100"
    return_code, sdterr, stdout = run_command(command)
    assert return_code == 0
    expected_files = {"database.dmnd", "database.fasta", "database.sqlite3"}
    confirm_files_present(expected_files, tmpdir)


def test_local_gne(tmpdir):
    command = \
        f"cblaster -d  --testing gne {TEST_FILE_DIR}{os.sep}test_session_local_gbk_{OS_NAME}.json --max_gap 200000" \
        f' --samples 25 --scale log -o {str(tmpdir.join("summary.txt"))} -hh -d "\t" -e 2 ' \
        f'-p {str(tmpdir.join("plot.html"))}'
    return_code, sdterr, stdout = run_command(command)
    assert return_code == 0
    expected_files = {"summary.txt", "plot.html"}
    confirm_files_present(expected_files, tmpdir)
    compare_file_pairs([["summary.txt",  "gne_local_summary.txt"]], tmpdir)


def test_remote_gne(tmpdir):
    command = \
        f"cblaster -d  --testing gne {TEST_FILE_DIR}{os.sep}test_session_remote_fa_{OS_NAME}.json --max_gap 250000" \
        f' --samples 25 --scale linear -o {str(tmpdir.join("summary.txt"))} -hh -d "\t" -e 2 ' \
        f'-p {str(tmpdir.join("plot.html"))}'
    return_code, sdterr, stdout = run_command(command)
    assert return_code == 0
    expected_files = {"summary.txt", "plot.html"}
    confirm_files_present(expected_files, tmpdir)
    compare_file_pairs([["summary.txt",  "gne_remote_summary.txt"]], tmpdir)


def test_local_extract(tmpdir):
    command = \
        f"cblaster -d --testing extract {TEST_FILE_DIR}{os.sep}test_session_local_gbk_{OS_NAME}.json "\
        f" -q AEK75493.1 AEK75491.1 AEK75490.1 -o {str(tmpdir.join('summary.txt'))} -sc JF752342.1:1-70000 -de _ -es" \
        f" -no"
    return_code, sdterr, stdout = run_command(command)
    assert return_code == 0
    expected_files = {"summary.txt"}
    confirm_files_present(expected_files, tmpdir)
    compare_file_pairs([["summary.txt",  "extract_local_summary.txt"]], tmpdir)


def test_remote_extract(tmpdir):
    command = \
        f"cblaster -d --testing extract {TEST_FILE_DIR}{os.sep}test_session_remote_fa_{OS_NAME}.json " \
        f" -o {str(tmpdir.join('summary.txt'))} -de _ -es -no -or Verrucosispora.*"
    return_code, sdterr, stdout = run_command(command)
    assert return_code == 0
    expected_files = {"summary.txt"}
    confirm_files_present(expected_files, tmpdir)
    compare_file_pairs([["summary.txt",  "extract_remote_summary.txt"]], tmpdir)


def test_local_extract_clusters(tmpdir):
    command = \
        f"cblaster -d --testing extract_clusters {TEST_FILE_DIR}{os.sep}test_session_local_gbk_{OS_NAME}.json -c 2 -o" \
        f" {tmpdir} -pf test_ -f genbank -mc 5"
    return_code, sdterr, stdout = run_command(command)
    assert return_code == 0
    expected_files = {"test_cluster2.gbk"}
    confirm_files_present(expected_files, tmpdir)
    compare_file_pairs([["test_cluster2.gbk",  "extract_clusters_cluster2_local_genbank.gbk"]], tmpdir)


def test_remote_extract_clusters(tmpdir):
    command = \
        f"cblaster -d --testing extract_clusters {TEST_FILE_DIR}{os.sep}test_session_remote_fa_{OS_NAME}.json " \
        f"-o {tmpdir} -or Verrucosispora.* -sc KF826676.1:5000-30000 -pf test_ -f genbank -st 20 -f bigscape"
    return_code, sdterr, stdout = run_command(command)
    assert return_code == 0
    expected_files = {"test_cluster11.gbk"}
    confirm_files_present(expected_files, tmpdir)
    compare_file_pairs([["test_cluster11.gbk",  "extract_clusters_cluster11_remote_bigscape.gbk"]], tmpdir)


def test_local_plot_clusters(tmpdir):
    command = \
        f"cblaster -d --testing plot_clusters {TEST_FILE_DIR}{os.sep}test_session_local_gbk_{OS_NAME}.json -c 2 -o" \
        f" {str(tmpdir.join('plot.html'))} -mc 5"
    return_code, sdterr, stdout = run_command(command)
    assert return_code == 0


def test_remote_plot_clusters(tmpdir):
    command = \
        f"cblaster -d --testing plot_clusters {TEST_FILE_DIR}{os.sep}test_session_remote_fa_{OS_NAME}.json -o" \
        f" {str(tmpdir.join('plot.html'))} -mc 5 -or Verrucosispora.* -sc KF826676.1:5000-30000 -st 20"
    return_code, sdterr, stdout = run_command(command)
    assert return_code == 0


if __name__ == '__main__':
    pytest.main(args=['-s', os.path.abspath(__file__)])

