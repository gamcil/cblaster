"""
Test suite for database.py
"""

import subprocess

from pathlib import Path

import pytest

from cblaster import database

TEST_DIR = Path(__file__).resolve().parent


def test_diamond_makedb(mocker):
    mocker.patch("cblaster.helpers.get_program_path", return_value="test_path")
    mocker.patch("subprocess.run")

    database.diamond_makedb("fasta", "name")
    subprocess.run.assert_called_once_with(
        ["test_path", "makedb", "--in", "fasta", "--db", "name"],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
