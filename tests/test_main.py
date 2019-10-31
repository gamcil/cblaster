#!/usr/bin/env python3

"""
Test suite for main.py
"""

import pytest
import pytest_mock

from clusterblaster import local, context, main, classes, remote


class MockOrganism(classes.Organism):
    def count_hit_clusters(self):
        return 1

    def summary(self):
        return "mocked"


def test_summarise(capsys, tmp_path):
    def return_text(x):
        return "test"

    organisms = [
        MockOrganism(name="test", strain="123"),
        MockOrganism(name="test2", strain="456"),
    ]

    # Test output file handle
    file = tmp_path / "test.txt"
    with file.open("w") as handle:
        main.summarise(organisms, output=handle)
    assert file.read_text() == "mocked\n\n\nmocked"

    # Test stdout
    main.summarise(organisms)
    captured = capsys.readouterr()
    assert captured.out == "mocked\n\n\nmocked"


def test_clusterblaster(mocker, tmp_path):

    mocker.patch("clusterblaster.local.search")
    mocker.patch("clusterblaster.remote.search")
    mocker.patch("clusterblaster.context.search")
    mocker.patch("clusterblaster.main.summarise")

    file = tmp_path / "test.txt"

    with file.open("w") as handle:
        main.clusterblaster(query_ids=["seq1"], mode="local", output=handle)

    main.clusterblaster(query_ids=["seq1"], mode="remote")

    local.search.assert_called_once()
    remote.search.assert_called_once()

    context.search.call_count = 2
    main.summarise.call_count = 2


def test_get_arguments_remote_defaults():
    assert vars(main.get_arguments(["-qf", "test"])) == {
        "output": None,
        "debug": False,
        "query_ids": None,
        "query_file": "test",
        "mode": "remote",
        "database": "nr",
        "entrez_query": None,
        "rid": None,
        "gap": 20000,
        "conserve": 3,
        "max_evalue": 0.01,
        "min_identity": 30,
        "min_coverage": 50,
    }


def test_get_arguments_remote_invalid_db():
    with pytest.raises(ValueError):
        main.get_arguments(["-qf", "test", "-m", "remote", "-db", "fake"])


def test_get_arguments_local_invalid_arg():
    with pytest.raises(ValueError):
        main.get_arguments(["-qf", "test", "-m", "local", "-eq", "entrez"])
        main.get_arguments(["-qf", "test", "-m", "local", "--rid", "rid"])
