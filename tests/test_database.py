"""
Test suite for database.py
"""

import subprocess

import pytest

from clusterblaster import database


def test_get_genbank_paths(tmp_path):
    d = tmp_path / "folder"
    d.mkdir()

    for ext in [".gb", ".gbk", ".genbank", ".fake", ""]:
        _p = d / f"test{ext}"
        _p.write_text("test")

    paths = database.get_genbank_paths(d)

    assert paths == [d / "test.genbank", d / "test.gb", d / "test.gbk"]


def test_diamond_makedb(mocker):
    mocker.patch("clusterblaster.helpers.get_program_path", return_value="test_path")
    mocker.patch("subprocess.run")

    database.diamond_makedb("fasta", "name")
    subprocess.run.assert_called_once_with(
        ["test_path", "makedb", "--in", "fasta", "--db", "name"],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
