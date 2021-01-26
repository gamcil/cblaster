
"""Simple testing functions simply making sure that certain commands can run with certain values

Warning: File comparissons and summary output rely heavily on diamond. These tests are made with diamond v 2.0.6
earlier or later versions will probably produce slightly different results and the provided databases will only
work with this exact version of diamond. That is why these versions are included"""

import os
import shutil
import time
from pathlib import Path
from tempfile import mkdtemp
import platform
import sys

COUNT = 0
TOTAL = 6
RUN = {"local": False, "remote": True}


def run_search_command(command, actual_expected=None):
    global COUNT
    COUNT += 1
    print(f"Running command {COUNT}/{TOTAL}: '{command}'")
    start_time = time.time()
    return_code = os.system(command)
    end_time = time.time()
    if return_code != 0:
        raise SystemExit(f"Command terminated with non zero exit status.")
    print(f"Command finished in {end_time - start_time:.3f} seconds.")
    if actual_expected is None:
        print()
        return
    for actual_file_path, expected_file_path in actual_expected:
        try:
            compare_files(actual_file_path, expected_file_path)
        except FileNotFoundError:
            raise AssertionError(f"File {actual_file_path} has not been created during the execution of: '{command}'.")
    print(f"All output files are as expected.\n")


def compare_files(actual_file_path, expected_file_path):
    with open(out_dir + os.sep + actual_file_path, "r") as actual_file, \
            open(comparrison_file_dir + os.sep + expected_file_path, "r") as expected_file:
        for line_index, actual_expected in enumerate(zip(actual_file, expected_file)):
            actual_line, expected_line = actual_expected
            if actual_line != expected_line:
                raise AssertionError(f"File {actual_file_path} and {expected_file_path} do not match at line "
                                     f"{line_index + 1}")


if __name__ == '__main__':
    os_name = platform.system().lower()
    if os_name == "darwin":
        raise SystemExit("Test for Mac OS are not yet configured.")

    out_dir = mkdtemp()
    # make sure the paths point the right way
    current_dir = str(Path(__file__).resolve().parent)
    test_file_dir = current_dir + os.sep + "test_files"
    comparrison_file_dir = current_dir + os.sep + "comparisson_files"
    os.chdir(current_dir)

    # add diamond to path
    old_path = os.environ["PATH"]
    os.environ["PATH"] = f"{str(Path('./diamond_files').resolve().absolute())}{os.pathsep}{old_path}"

    try:
        # LOCAL TESTS
        if RUN["local"]:
            # test gbk query in local mode with all options enabled
            command1 = f"cblaster -d search -m local -qf {test_file_dir}{os.sep}test_query.gb -o " \
                       f"{out_dir}{os.sep}summary.txt -db {test_file_dir}{os.sep}test_database_{os_name}.dmnd -ohh" \
                       f" -ode , -odc 2 -osc -b {out_dir}{os.sep}binary.txt -bhh -bde _ -bdc 2 -bkey sum -bat coverage " \
                       f" --blast_file {out_dir}{os.sep}blast.txt --ipg_file {out_dir}{os.sep}ipgs.txt " \
                       f"-g 25000 -u 2 -mh 3 -r AEK75493.1 -me 0.01 -mi 30 -mc 50 -s {out_dir}{os.sep}session.json"
            actual_vs_expected_files = [["summary.txt", "summary_local_gbk.txt"], ["binary.txt", "binary_local_gbk.txt"]]
            run_search_command(command1, actual_vs_expected_files)
            os.remove(f"{out_dir}{os.sep}session.json")

            # test embl query in local mode with all options enabled
            command2 = f"cblaster -d search -m local -qf {test_file_dir}{os.sep}test_query.embl -o " \
                       f"{out_dir}{os.sep}summary.txt -db {test_file_dir}{os.sep}test_database_{os_name}.dmnd -ohh" \
                       f" -ode , -odc 2 -osc -b {out_dir}{os.sep}binary.txt -bhh -bde _ -bdc 2 -bkey sum -bat coverage " \
                       f" --blast_file {out_dir}{os.sep}blast.txt --ipg_file {out_dir}{os.sep}ipgs.txt " \
                       f"-g 25000 -u 2 -mh 3 -r AEK75493.1 -me 0.01 -mi 30 -mc 50 -s {out_dir}{os.sep}session.json"
            run_search_command(command2)
            os.remove(f"{out_dir}{os.sep}session.json")

            # test fasta query in local mode with all options enabled
            command3 = f"cblaster -d search -m local -qf {test_file_dir}{os.sep}test_query.fa -o " \
                       f"{out_dir}{os.sep}summary.txt -db {test_file_dir}{os.sep}test_database_{os_name}.dmnd -ohh" \
                       f" -ode , -odc 2 -osc -b {out_dir}{os.sep}binary.txt -bhh -bde _ -bdc 2 -bkey sum -bat coverage " \
                       f" --blast_file {out_dir}{os.sep}blast.txt --ipg_file {out_dir}{os.sep}ipgs.txt " \
                       f"-g 25000 -u 2 -mh 3 -r AEK75493.1 -me 0.01 -mi 30 -mc 50 -s {out_dir}{os.sep}session.json"
            run_search_command(command3)
            os.remove(f"{out_dir}{os.sep}session.json")

            # test query identifiers in local mode
            command4 = f"cblaster -d search -m local -qi AEK75490.1 AEK75490.1 AEK75500.1 AEK75516.1 AEK75516.1" \
                       f" AEK75502.1 -o {out_dir}{os.sep}summary.txt -db " \
                       f"{test_file_dir}{os.sep}test_database_{os_name}.dmnd -ohh -ode , -odc 2 -osc -b" \
                       f" {out_dir}{os.sep}binary.txt -bhh -bde _ -bdc 2 -bkey sum -bat coverage " \
                       f" --blast_file {out_dir}{os.sep}blast.txt --ipg_file {out_dir}{os.sep}ipgs.txt " \
                       f"-g 25000 -u 2 -mh 3 -me 0.01 -mi 30 -mc 50 -s {out_dir}{os.sep}session.json"
            run_search_command(command4)
            os.remove(f"{out_dir}{os.sep}session.json")

            # test local session with all options enabled
            command5 = f"cblaster -d search -m local -qf {test_file_dir}{os.sep}test_query.gb " \
                       f"-s {test_file_dir}{os.sep}test_session_local_embl.json " \
                       f"{test_file_dir}{os.sep}test_session_local_gbk.json -db" \
                       f" {test_file_dir}{os.sep}test_database_{os_name}.dmnd -o {out_dir}{os.sep}summary.txt"
            actual_vs_expected_files = [["summary.txt", "summary_local_gbk_embl_combined.txt"]]
            run_search_command(command5, actual_vs_expected_files)

            # test recompute
            shutil.copyfile(f"{test_file_dir}{os.sep}test_session_local_gbk.json",
                            f"{out_dir}{os.sep}test_session_local_gbk_copy.json")
            command6 = f"cblaster -d search -m local -qf {test_file_dir}{os.sep}test_query.gb " \
                       f"-s {out_dir}{os.sep}test_session_local_gbk_copy.json --recompute" \
                       f" -db {test_file_dir}{os.sep}test_database_{os_name}.dmnd -o {out_dir}{os.sep}summary.txt " \
                       f"-g 50000 -u 5 -mh 3 -me 0.01 -mi 30 -mc 50 -b {out_dir}{os.sep}binary.txt"
            actual_vs_expected_files = [["summary.txt", "summary_local_gbk_recompute.txt"],
                                        ["binary.txt", "binary_local_gbk_recompute.txt"]]
            run_search_command(command6, actual_vs_expected_files)

    # REMOTE TESTS --> all remote tests are conducted with a session file because the runtime of remote runs can
    # vary greatly
        if RUN["remote"]:
            # load a remote session and make summary and binary files
            command1 = f"cblaster -d search -m remote -qf {test_file_dir}{os.sep}test_query.fa -s " \
                       f"{test_file_dir}{os.sep}test_session_remote_fa.json" \
                       f" -o {out_dir}{os.sep}summary.txt -b {out_dir}{os.sep}binary.txt"
            actual_vs_expected_files = [["summary.txt", "summary_remote_fa.txt"],
                                        ["binary.txt", "binary_remote_fa.txt"]]
            run_search_command(command1, actual_vs_expected_files)

            # recompute a remote session
            shutil.copyfile(f"{test_file_dir}{os.sep}test_session_remote_fa.json",
                            f"{out_dir}{os.sep}test_session_remote_fa_copy.json")
            command2 = f"cblaster -d search -m remote -qf {test_file_dir}{os.sep}test_query.fa -s " \
                       f"{out_dir}{os.sep}test_session_remote_fa_copy.json" \
                       f" -o {out_dir}{os.sep}summary.txt -b {out_dir}{os.sep}binary.txt --recompute -g 50000 -u 7" \
                       f" -mh 3 -me 0.01 -mi 20 -mc 60 --sort_clusters"
            actual_vs_expected_files = [["summary.txt", "summary_remote_fa_recompute.txt"],
                                        ["binary.txt", "binary_remote_fa_recompute.txt"]]
            run_search_command(command2, actual_vs_expected_files)

    # make sure to always remove the dir even on error
    finally:
        shutil.rmtree(out_dir)
