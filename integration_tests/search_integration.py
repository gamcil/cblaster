
"""Simple testing functions simply making sure that certain commands can run with certain values"""

import os
import shutil
import time
from pathlib import Path
from tempfile import mkdtemp

COUNT = 0
TOTAL = 5


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
        for actual_line, expected_line in zip(actual_file, expected_file):
            if actual_line != expected_line:
                raise AssertionError(f"Expected {expected_line} but got {actual_line}. File {expected_file_path}"
                                     f" and {actual_file_path} do not match.")


if __name__ == '__main__':
    out_dir = mkdtemp()
    # make sure the paths point the right way
    current_dir = str(Path(__file__).resolve().parent)
    test_file_dir = current_dir + os.sep + "test_files"
    comparrison_file_dir = current_dir + os.sep + "comparisson_files"
    os.chdir(current_dir)

    # LOCAL TESTS
    try:
        # test gbk query in local mode with all options enabled
        command1 = f"cblaster -d search -m local -qf {test_file_dir}{os.sep}test_query.gb -o " \
                   f"{out_dir}{os.sep}summary.txt -db {test_file_dir}{os.sep}test_database.dmnd -ohh -ode , -odc 2 " \
                   f"-osc -b {out_dir}{os.sep}binary.txt -bhh -bde _ -bdc 2 -bkey sum -bat coverage " \
                   f" --blast_file {out_dir}{os.sep}blast.txt --ipg_file {out_dir}{os.sep}ipgs.txt " \
                   f"-g 25000 -u 2 -mh 3 -r AEK75493.1 -me 0.01 -mi 30 -mc 50 -s {out_dir}{os.sep}session.json"
        actual_vs_expected_files = [["summary.txt", "summary_local_gbk.txt"], ["binary.txt", "binary_local_gbk.txt"],
                                    ["blast.txt", "blast_local_gbk.txt"], ["session.json", "session_local_gbk.json"]]
        run_search_command(command1, actual_vs_expected_files)
        os.remove(f"{out_dir}{os.sep}session.json")

        # test embl query in local mode with all options enabled
        command2 = f"cblaster -d search -m local -qf {test_file_dir}{os.sep}test_query.embl -o " \
                   f"{out_dir}{os.sep}summary.txt -db {test_file_dir}{os.sep}test_database.dmnd -ohh -ode , -odc 2 " \
                   f"-osc -b {out_dir}{os.sep}binary.txt -bhh -bde _ -bdc 2 -bkey sum -bat coverage " \
                   f" --blast_file {out_dir}{os.sep}blast.txt --ipg_file {out_dir}{os.sep}ipgs.txt " \
                   f"-g 25000 -u 2 -mh 3 -r AEK75493.1 -me 0.01 -mi 30 -mc 50 -s {out_dir}{os.sep}session.json"
        run_search_command(command2)
        os.remove(f"{out_dir}{os.sep}session.json")

        # test fasta query in local mode with all options enabled
        command3 = f"cblaster -d search -m local -qf {test_file_dir}{os.sep}test_query.fa -o " \
                   f"{out_dir}{os.sep}summary.txt -db {test_file_dir}{os.sep}test_database.dmnd -ohh -ode , -odc 2 " \
                   f"-osc -b {out_dir}{os.sep}binary.txt -bhh -bde _ -bdc 2 -bkey sum -bat coverage " \
                   f" --blast_file {out_dir}{os.sep}blast.txt --ipg_file {out_dir}{os.sep}ipgs.txt " \
                   f"-g 25000 -u 2 -mh 3 -r AEK75493.1 -me 0.01 -mi 30 -mc 50 -s {out_dir}{os.sep}session.json"
        run_search_command(command3)
        os.remove(f"{out_dir}{os.sep}session.json")

        # test query identifiers in local mode
        command4 = f"cblaster -d search -m local -qi AEK75490.1 AEK75490.1 AEK75500.1 AEK75516.1 AEK75516.1" \
                   f" AEK75502.1 -o {out_dir}{os.sep}summary.txt -db {test_file_dir}{os.sep}test_database.dmnd -ohh " \
                   f"-ode , -odc 2 -osc -b {out_dir}{os.sep}binary.txt -bhh -bde _ -bdc 2 -bkey sum -bat coverage " \
                   f" --blast_file {out_dir}{os.sep}blast.txt --ipg_file {out_dir}{os.sep}ipgs.txt " \
                   f"-g 25000 -u 2 -mh 3 -me 0.01 -mi 30 -mc 50 -s {out_dir}{os.sep}session.json"
        run_search_command(command4)
        os.remove(f"{out_dir}{os.sep}session.json")

        # test local session with all options enabled
        command5 = f"cblaster -d search -m local -qf {test_file_dir}{os.sep}test_query.gb " \
                   f"-s {test_file_dir}{os.sep}test_session_embl.json {test_file_dir}{os.sep}test_session_gbk.json " \
                   f"-db {test_file_dir}{os.sep}test_database.dmnd -o {out_dir}{os.sep}summary.txt"
        actual_vs_expected_files = [["summary.txt", "summary_local_session_combined.txt"]]
        run_search_command(command5, actual_vs_expected_files)

    # make sure to always remove the dir even on error
    finally:
        shutil.rmtree(out_dir)
