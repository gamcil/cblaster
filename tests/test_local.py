#!/usr/bin/env python3

"""
Test suite for local.py
"""


import subprocess
import pytest

from pathlib import Path

from cblaster import local, helpers


TEST_DIR = Path(__file__).resolve().parent


def test_parse():
    # length of QBE85648 == 179
    result = [
        # qid sid pid len mismatch gapopen qstart qend sstart ssend evalue bitscore
        "QBE85648.1\tHIT1\t100.000\t100.000\t1.38e-127\t365",
        "QBE85648.1\tHIT2\t20.000\t100.000\t1.38e-127\t365",
        "QBE85648.1\tHIT3\t100.000\t20.000\t1.38e-127\t365",
        "QBE85648.1\tHIT4\t100.000\t100.000\t0.011\t365",
    ]

    hits = local.parse(result)

    # Default thresholds are 30% identity, 50% coverage, 0.01 evalue
    # so only the first hit should be saved
    assert len(hits) == 1
    assert hits[0].query == "QBE85648.1"
    assert hits[0].subject == "HIT1"
    assert hits[0].identity == 100.0
    assert hits[0].coverage == 100.0
    assert hits[0].bitscore == 365.0
    assert hits[0].evalue == 1.38e-127


def test_parse_no_hits():
    with pytest.raises(SystemExit):
        local.parse([])


def test_diamond(monkeypatch):
    def mock_path(aliases):
        return "diamond"

    def mock_run(command, **kwargs):
        assert command == [
            "diamond",
            "blastp",
            "--query",
            "fasta",
            "--db",
            "database",
            "--id",
            "30",
            "--evalue",
            "0.01",
            "--outfmt",
            "6",
            "qseqid",
            "sseqid",
            "pident",
            "qcovhsp",
            "evalue",
            "bitscore",
            "--threads",
            "1",
            "--query-cover",
            "50",
            "--max-hsps",
            "1",
        ]
        return subprocess.CompletedProcess(
            args=command, stdout=b"line1\nline2\nline3", returncode=1
        )

    monkeypatch.setattr(helpers, "get_program_path", mock_path)
    monkeypatch.setattr(subprocess, "run", mock_run)

    assert local.diamond("fasta", "database") == ["line1", "line2", "line3"]


def test_search_ids(monkeypatch):
    def mock_efetch(ids):
        return {"SEQ1": "ABCDEF", "SEQ2": "ABCDEF", "SEQ3": "ABCDEF"}

    def mock_search(file, db, **kwargs):
        with open(file) as handle:
            text = handle.read()
        assert text == ">SEQ1\nABCDEF\n>SEQ2\nABCDEF\n>SEQ3\nABCDEF\n"

    monkeypatch.setattr(helpers, "efetch_sequences", mock_efetch)
    monkeypatch.setattr(local, "_search_file", mock_search)

    local._search_ids(["SEQ1", "SEQ2", "SEQ3"], "database")


def test_search_no_input():
    with pytest.raises(ValueError):
        local.search(database="database")


def test_search_query_ids(mocker):
    mocker.patch("cblaster.local._search_ids")
    local.search("database", query_ids=["test"])
    local._search_ids.assert_called_once_with(["test"], "database")


def test_search_query_files(mocker):
    mocker.patch("cblaster.local._search_file")
    local.search("database", query_file="test")
    local._search_file.assert_called_once_with("test", "database")
