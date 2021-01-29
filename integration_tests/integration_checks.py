"""Simple testing functions simply making sure that commands can run with certain values and some basic output checks

Warning: File comparissons and summary output rely heavily on diamond. These tests are made with diamond v 2.0.6
earlier or later versions will probably produce slightly different results and the provided databases will only
work with this exact version of diamond. That is why these versions are included"""

import os
import subprocess
import sys
import time
from pathlib import Path
from tempfile import mkdtemp
import platform

OUT_DIR = None
CURRENT_DIR = None
TEST_FILE_DIR = None
COMPARISSON_FILE_DIR = None
OS_NAME = platform.system().lower()

POSSIBLE_FLAGS = ["silent"]


class CommandTest:
    def __init__(self, command, description, compare_files=None):
        self.command = command
        self.run_time = 0
        self.description = description
        self.return_code = None
        self.__actual_expected = compare_files if compare_files is not None else []

    def run(self, silent=False):
        global OUT_DIR
        try:
            self.__run_command(silent)
        finally:
            # delete all files but not the folder
            for file_name in os.listdir(OUT_DIR):
                try:
                    file_path = os.path.join(OUT_DIR, file_name)
                    os.unlink(file_path)
                except Exception as e:
                    print(f"Failed to delete {file_path}. Reason:{e}")

    def __run_command(self, silent=False):
        print(f"Running command: '{self.description}'")
        start_time = time.time()
        kwargs = dict(stdout=subprocess.DEVNULL)
        if silent:
            kwargs["stderr"] = subprocess.DEVNULL
        popen = subprocess.Popen(self.command, shell=True, **kwargs)
        stdout, stderer = popen.communicate()
        end_time = time.time()
        self.return_code = popen.returncode
        if self.return_code != 0:
            raise RuntimeError(f"WARNING: command {self.description} failed.")
        self.run_time = end_time - start_time
        print(f"Command finished in {self.run_time:.3f} seconds.")
        if not self.__actual_expected:
            print()
            return
        if self.return_code != 0:
            print()
            return
        for actual_file_path, expected_file_path in self.__actual_expected:
            try:
                self.compare_files(actual_file_path, expected_file_path)
            except FileNotFoundError:
                raise AssertionError(
                    f"File {actual_file_path} has not been created during the execution of: '{self.description}'.")
        print(f"All output files are as expected.\n")

    def compare_files(self, actual_file_path, expected_file_path):
        global OUT_DIR, COMPARISSON_FILE_DIR
        with open(OUT_DIR + os.sep + actual_file_path, "r") as actual_file, \
                open(COMPARISSON_FILE_DIR + os.sep + expected_file_path, "r") as expected_file:
            for line_index, actual_expected in enumerate(zip(actual_file, expected_file)):
                actual_line, expected_line = actual_expected
                if actual_line != expected_line:
                    raise AssertionError(f"File {actual_file_path} and {expected_file_path} do not match at line "
                                         f"{line_index + 1}")


def setup():
    global OS_NAME
    if OS_NAME == "darwin":
        raise SystemExit("Test for Mac OS are not yet configured.")

    global OUT_DIR, TEST_FILE_DIR, COMPARISSON_FILE_DIR, CURRENT_DIR
    OUT_DIR = mkdtemp()
    # make sure the paths point the right way
    CURRENT_DIR = str(Path(__file__).resolve().parent)
    TEST_FILE_DIR = CURRENT_DIR + os.sep + "test_files"
    COMPARISSON_FILE_DIR = CURRENT_DIR + os.sep + "comparisson_files"
    os.chdir(CURRENT_DIR)

    # add the local versions of diamond to path
    old_path = os.environ["PATH"]
    os.environ["PATH"] = f"{str(Path('./diamond_files').resolve().absolute())}{os.pathsep}{old_path}"


def test_commands(arguments):
    flags, command_names = filter_flags(arguments)
    commands = []
    command_names = command_names if command_names else ["search", "makedb", "gne", "extract", "extract_clusters"]
    for name in command_names:
        name = name.lower()
        if name == "search":
            commands.extend(search_commands())
        elif name == "search_local":
            commands.extend(search_local_commands())
        elif name == "search_remote":
            commands.extend(search_remote_commands())
        elif name == "makedb":
            commands.extend(makedb_commands())
        elif name == "gne":
            commands.extend(gne_commands())
        elif name == "extract":
            commands.extend(extract_commands())
        elif name == "extract_clusters":
            commands.extend(extract_clusters_commands())
        else:
            raise ValueError(f"No command named {name}.")
    for command in commands:
        command.run(flags["silent"])
    print("All tests passed!")


def filter_flags(arguments):
    flags = {name: False for name in POSSIBLE_FLAGS}
    command_names = []
    for item in arguments:
        if item in POSSIBLE_FLAGS:
            flags[item] = True
        else:
            command_names.append(item)
    return flags, command_names


def search_commands():
    global OUT_DIR, TEST_FILE_DIR
    return search_local_commands() + search_remote_commands()


def search_local_commands():
    global OUT_DIR, TEST_FILE_DIR

    commands = [
        # test gbk query in local mode with all options enabled
        CommandTest(
            f"cblaster -d search -m local -qf {TEST_FILE_DIR}{os.sep}test_query.gb -o "
            f"{OUT_DIR}{os.sep}summary.txt -db {TEST_FILE_DIR}{os.sep}test_database_{OS_NAME}.dmnd -ohh"
            f" -ode , -odc 2 -osc -b {OUT_DIR}{os.sep}binary.txt -bhh -bde _ -bdc 2 -bkey sum -bat coverage "
            f" --blast_file {OUT_DIR}{os.sep}blast.txt --ipg_file {OUT_DIR}{os.sep}ipgs.txt "
            f"-g 25000 -u 2 -mh 3 -r AEK75493.1 -me 0.01 -mi 30 -mc 50 -s {OUT_DIR}{os.sep}session.json",
            "test gbk query local",
            [["summary.txt", "summary_local_gbk.txt"], ["binary.txt", "binary_local_gbk.txt"]]
        ),
        # test embl query in local mode with all options enabled
        CommandTest(
            f"cblaster -d search -m local -qf {TEST_FILE_DIR}{os.sep}test_query.embl -o "
            f"{OUT_DIR}{os.sep}summary.txt -db {TEST_FILE_DIR}{os.sep}test_database_{OS_NAME}.dmnd -ohh"
            f" -ode , -odc 2 -osc -b {OUT_DIR}{os.sep}binary.txt -bhh -bde _ -bdc 2 -bkey sum -bat coverage "
            f" --blast_file {OUT_DIR}{os.sep}blast.txt --ipg_file {OUT_DIR}{os.sep}ipgs.txt "
            f"-g 25000 -u 2 -mh 3 -r AEK75493.1 -me 0.01 -mi 30 -mc 50 -s {OUT_DIR}{os.sep}session.json",
            "test embl query local"
        ),
        # test fasta query in local mode with all options enabled
        CommandTest(
            f"cblaster -d search -m local -qf {TEST_FILE_DIR}{os.sep}test_query.fa -o "
            f"{OUT_DIR}{os.sep}summary.txt -db {TEST_FILE_DIR}{os.sep}test_database_{OS_NAME}.dmnd -ohh"
            f" -ode , -odc 2 -osc -b {OUT_DIR}{os.sep}binary.txt -bhh -bde _ -bdc 2 -bkey sum -bat coverage "
            f" --blast_file {OUT_DIR}{os.sep}blast.txt --ipg_file {OUT_DIR}{os.sep}ipgs.txt "
            f"-g 25000 -u 2 -mh 3 -r AEK75493.1 -me 0.01 -mi 30 -mc 50 -s {OUT_DIR}{os.sep}session.json",
            "test fasta query local"
        ),
        # test query identifiers in local mode
        CommandTest(
            f"cblaster -d search -m local -qi AEK75490.1 AEK75490.1 AEK75500.1 AEK75516.1 AEK75516.1"
            f" AEK75502.1 -o {OUT_DIR}{os.sep}summary.txt -db "
            f"{TEST_FILE_DIR}{os.sep}test_database_{OS_NAME}.dmnd -ohh -ode , -odc 2 -osc -b"
            f" {OUT_DIR}{os.sep}binary.txt -bhh -bde _ -bdc 2 -bkey sum -bat coverage "
            f" --blast_file {OUT_DIR}{os.sep}blast.txt --ipg_file {OUT_DIR}{os.sep}ipgs.txt "
            f"-g 25000 -u 2 -mh 3 -me 0.01 -mi 30 -mc 50 -s {OUT_DIR}{os.sep}session.json",
            "test query identifiers local"
        ),
        # test local session with all options enabled
        CommandTest(
            f"cblaster -d search -m local -qf {TEST_FILE_DIR}{os.sep}test_query.gb "
            f"-s {TEST_FILE_DIR}{os.sep}test_session_local_embl_{OS_NAME}.json "
            f"{TEST_FILE_DIR}{os.sep}test_session_local_gbk_{OS_NAME}.json -db"
            f" {TEST_FILE_DIR}{os.sep}test_database_{OS_NAME}.dmnd -o {OUT_DIR}{os.sep}summary.txt",
            "test local session",
            [["summary.txt", "summary_local_gbk_embl_combined.txt"]]
        ),
        # test recompute a local session
        CommandTest(
            f"cblaster -d search -m local -qf {TEST_FILE_DIR}{os.sep}test_query.gb "
            f"-s {OUT_DIR}{os.sep}test_session_local_gbk_copy.json --recompute"
            f" -db {TEST_FILE_DIR}{os.sep}test_database_{OS_NAME}.dmnd -o {OUT_DIR}{os.sep}summary.txt "
            f"-g 50000 -u 5 -mh 3 -me 0.01 -mi 30 -mc 50 -b {OUT_DIR}{os.sep}binary.txt",
            "test recompute local",
            [["summary.txt", "summary_local_gbk_recompute.txt"], ["binary.txt", "binary_local_gbk_recompute.txt"]]
        )]
    return commands


def search_remote_commands():
    global OUT_DIR, TEST_FILE_DIR

    # load a remote session and make summary and binary files
    commands = [
        CommandTest(
            f"cblaster -d search -m remote -qf {TEST_FILE_DIR}{os.sep}test_query.fa -s "
            f"{TEST_FILE_DIR}{os.sep}test_session_remote_fa.json"
            f" -o {OUT_DIR}{os.sep}summary.txt -b {OUT_DIR}{os.sep}binary.txt",
            "load remote session",
            [["summary.txt", "summary_remote_fa.txt"], ["binary.txt", "binary_remote_fa.txt"]]
        ),
        # recompute a remote session
        CommandTest(
            f"cblaster -d search -m remote -qf {TEST_FILE_DIR}{os.sep}test_query.fa -s "
            f"{TEST_FILE_DIR}{os.sep}test_session_remote_fa.json"
            f" -o {OUT_DIR}{os.sep}summary.txt -b {OUT_DIR}{os.sep}binary.txt --recompute -g 50000 -u 7"
            f" -mh 3 -me 0.01 -mi 20 -mc 60 --sort_clusters",
            "recompute remote session",
            [["summary.txt", "summary_remote_fa_recompute.txt"], ["binary.txt", "binary_remote_fa_recompute.txt"]]
        )
    ]
    return commands


def makedb_commands():
    global OUT_DIR, TEST_FILE_DIR
    commands = [
        # using embl and gbk files
        CommandTest(
            f"cblaster -d makedb {TEST_FILE_DIR}{os.sep}test_query.gb"
            f" {TEST_FILE_DIR}{os.sep}test_query.embl -n {OUT_DIR}{os.sep}database -b 3 -c 100",
            "makedb gbk, embl files",
            [["database.fasta", "makedb_database_gbk_embl.fasta"]]
        ),
        # using gff files
        CommandTest(
            f"cblaster -d makedb {TEST_FILE_DIR}{os.sep}test_gff_v_maris.fna"
            f" {TEST_FILE_DIR}{os.sep}test_gff_v_maris.gff -n {OUT_DIR}{os.sep}database -b 3 -c 100",
            "makedb gff files"
        )
    ]
    return commands


def gne_commands():
    global OUT_DIR, TEST_FILE_DIR
    commands = [
        CommandTest(
            f"cblaster -d gne {TEST_FILE_DIR}{os.sep}test_session_local_gbk_{OS_NAME}.json --max_gap 200000 --samples 25"
            f' --scale log -o {OUT_DIR}{os.sep}summary.txt -hh -d "\t" -e 2',
            "gne with local gbk session",
            [["summary.txt", "gne_local_summary.txt"]]
        ),
        CommandTest(
            f"cblaster -d gne {TEST_FILE_DIR}{os.sep}test_session_remote_fa.json --max_gap 250000 --samples 10"
            f' --scale linear -o {OUT_DIR}{os.sep}summary.txt -hh -d "\t" -e 3',
            "gne with remote fa session",
            [["summary.txt", "gne_remote_summary.txt"]]
        )
    ]
    return commands


def extract_commands():
    global OUT_DIR, TEST_FILE_DIR
    commands = [
        CommandTest(
            f"cblaster -d extract {TEST_FILE_DIR}{os.sep}test_session_local_gbk_{OS_NAME}.json "
            f"-or GCA_0002041.* -q AEK75517.1 AEK75496.1 AEK75490.1 -o"
            f" {OUT_DIR}{os.sep}summary.txt -sc CP002638.1:2562030-2619476 -de _ -es -no",
            "local session extraction",
            [["summary.txt", "extract_local_summary.txt"]]
        ),
        CommandTest(
            f"cblaster -d extract {TEST_FILE_DIR}{os.sep}test_session_remote_fa.json -o "
            f"{OUT_DIR}{os.sep}summary.txt -de : -es -or Verrucosispora.*",
            "remote session extraction",
            [["summary.txt", "extract_remote_summary.txt"]]
        )
    ]
    return commands


def extract_clusters_commands():
    global OUT_DIR, TEST_FILE_DIR
    commands = [
        CommandTest(
            f"cblaster -d extract_clusters {TEST_FILE_DIR}{os.sep}test_session_local_gbk_{OS_NAME}.json -c 1-2 4 -o {OUT_DIR}"
            f" -pf test_ -f genbank -mec 5",
            "local session cluster extraction",
            [["test_cluster1.gbk", "extract_clusters_cluster1_local_genbank.gbk"]]
        ),
        CommandTest(
            f"cblaster -d extract_clusters {TEST_FILE_DIR}{os.sep}test_session_remote_fa.json -o {OUT_DIR} -or"
            f" Verrucosispora.* -sc KF826676.1:5000-30000 -st 20 -f bigscape",
            "remote session cluster extraction",
            [["cluster11.gbk", "extract_clusters_cluster11_remote_bigscape.gbk"]]
        )
    ]
    return commands


if __name__ == '__main__':
    # run tests of all the requested commands by specifying the command name separated by spaces
    cmd_arguments = sys.argv[1:]
    print("Running setup")
    setup()
    print("\nRunning tests:\n")
    test_commands(cmd_arguments)
