
import pytest


# for testing purposes
import sys
sys.path.insert(0, "..")

from cblaster import parsers


def test_subcommands():
    # valid subcommands
    for subcommand in [["gui"],
                       ["makedb", "placeholder_genbanks", "placeholder_filename"],
                       ["gne", "placeholder_session"],
                       ["extract", "placeholder_session"],
                       ["search"]]:
        parsers.parse_args(subcommand)

    # an invalid subcommand
    with pytest.raises(SystemExit):
        parsers.parse_args(["test"])


def test_remote_db():
    # make sure that no error is raised with these database types
    for name in ["nr", "refseq_protein", "swissprot", "pdbaa"]:
        parsers.parse_args(["search", "-qf", "test.faa", "-m", "remote", "-db", name])
    # make sure an error is raised on invalid database
    with pytest.raises(SystemExit):
        parsers.parse_args(["search", "-qf", "test.faa", "-m", "remote", "-db", "test"])


def test__local_invalid_arg():
    with pytest.raises(SystemExit):
        parsers.parse_args(["search", "-qf", "test.faa", "-m", "local", "-eq", "entrez"])
        parsers.parse_args(["search", "-qf", "test.faa", "-m", "local", "--rid", "rid"])
