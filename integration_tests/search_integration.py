
"""Simple testing functions simply making sure that certain commands can run with certain values"""

import os
import shutil
import time
from pathlib import Path
from tempfile import mkdtemp

COUNT = 0
TOTAL = 3


def search_run_test(command):
    global COUNT
    COUNT += 1
    print(f"Running command {COUNT}/{TOTAL}: '{command}'")
    start_time = time.time()
    return_code = os.system(command)
    end_time = time.time()
    if return_code != 0:
        raise SystemExit(f"Test search command '{command}' terminated with non zero exit status")
    print(f"command finished in {end_time - start_time:.3f} seconds")


if __name__ == '__main__':
    out_dir = mkdtemp()
    # make sure the paths point the right way
    current_dir = str(Path(__file__).resolve().parent)
    test_file_dir = current_dir + os.sep + "test_files"
    os.chdir(current_dir)

    # test gbk query and as much options as possible
    command1 = f"cblaster -d search -m local -qf {test_file_dir}{os.sep}test_query.gb -o {out_dir}{os.sep}summary.txt" \
               f" -db {test_file_dir}{os.sep}test_database.dmnd -ohh -ode , -odc 2 -osc -b " \
               f"{out_dir}{os.sep}binary.txt -bhh -bde , -bdc 2 -bkey sum -bat coverage " \
               f" --blast_file {out_dir}{os.sep}blast_save.txt --ipg_file {out_dir}{os.sep}ipgs.txt " \
               f"-g 25000 -u 2 -mh 3 -r AEK75493.1 -me 0.01 -mi 30 -mc 50"
    search_run_test(command1)

    # test embl query and as much options as possible
    command2 = f"cblaster -d search -m local -qf {test_file_dir}{os.sep}test_query.embl -o {out_dir}{os.sep}summary.txt" \
               f" -db {test_file_dir}{os.sep}test_database.dmnd -ohh -ode , -odc 2 -osc -b " \
               f"{out_dir}{os.sep}binary.txt -bhh -bde , -bdc 2 -bkey sum -bat coverage " \
               f" --blast_file {out_dir}{os.sep}blast_save.txt --ipg_file {out_dir}{os.sep}ipgs.txt " \
               f"-g 25000 -u 2 -mh 3 -r AEK75493.1 -me 0.01 -mi 30 -mc 50"
    search_run_test(command2)

    # test fasta query and as much options as possible
    command2 = f"cblaster -d search -m local -qf {test_file_dir}{os.sep}test_query.fa -o {out_dir}{os.sep}summary.txt" \
               f" -db {test_file_dir}{os.sep}test_database.dmnd -ohh -ode , -odc 2 -osc -b " \
               f"{out_dir}{os.sep}binary.txt -bhh -bde , -bdc 2 -bkey sum -bat coverage " \
               f" --blast_file {out_dir}{os.sep}blast_save.txt --ipg_file {out_dir}{os.sep}ipgs.txt " \
               f"-g 25000 -u 2 -mh 3 -r AEK75493.1 -me 0.01 -mi 30 -mc 50"
    search_run_test(command2)

    shutil.rmtree(out_dir)
