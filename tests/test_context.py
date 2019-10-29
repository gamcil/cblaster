#!/usr/bin/env python3

"""
Test suite for context.py
"""

import pytest

from pathlib import Path

import requests
import requests_mock

from clusterblaster import classes, context


TEST_DIR = Path(__file__).resolve().parent


@pytest.fixture()
def hits():
    return [
        classes.Hit("q1", "s1", "1", "1", "0", "1", 0, 1000, "+"),
        classes.Hit("q2", "s2", "1", "1", "0", "1", 2000, 3000, "+"),
        classes.Hit("q3", "s3", "1", "1", "0", "1", 5000, 6000, "+"),
        classes.Hit("q4", "s4", "1", "1", "0", "1", 9000, 10000, "+"),
        classes.Hit("q5", "s5", "1", "1", "0", "1", 14000, 15000, "+"),
    ]


@pytest.mark.parametrize(
    "conserve, gap, results",
    [
        (0, 999, [[0], [1], [2], [3], [4]]),
        (0, 1000, [[0, 1], [2], [3], [4]]),
        (2, 2000, [[0, 1, 2]]),
        (4, 2000, []),
        (6, 2000, []),
    ],
)
def test_find_clusters(hits, conserve, gap, results):
    groups = context.find_clusters_in_hits(hits, conserve=conserve, gap=gap)

    assert len(groups) == len(results)

    for group, result in zip(groups, results):
        assert group == [hits[i] for i in result]


@pytest.mark.parametrize("conserve, gap", [(-1, 100), (1, -1)])
def test_find_clusters_negative_input(hits, conserve, gap):
    with pytest.raises(ValueError):
        context.find_clusters_in_hits(hits, conserve=conserve, gap=gap)


def test_efetch_IPGs(hits):
    with requests_mock.Mocker() as mock:
        mock.post("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?")

        context.efetch_IPGs([hit.subject for hit in hits])

        assert mock.request_history[0].url == (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
            "db=protein"
            "&rettype=ipg"
            "&retmode=text"
        )


def test_efetch_IPGs_error(hits):
    with requests_mock.Mocker() as mock, pytest.raises(requests.HTTPError):
        mock.post(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?",
            status_code=400,
        )
        context.efetch_IPGs([hit.subject for hit in hits])


def test_efetch_IPGs_output(hits, tmp_path):
    with requests_mock.Mocker() as mock:
        mock.post(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?", text="test"
        )

        test_out = tmp_path / "out.tsv"

        with test_out.open("w") as handle:
            context.efetch_IPGs([hit.subject for hit in hits], output_handle=handle)

        assert test_out.read_text() == "test"


def test_parse_IPG_table(hits):
    results = TEST_DIR / "ipg_results.txt"

    with open(results) as handle:
        parsed = context.parse_IPG_table(handle, hits)

    assert len(parsed) == 2

    one, two = parsed

    assert one.name == two.name == "Test organism"
    assert one.strain == "STRAIN 1"
    assert two.strain == "STRAIN 2"

    assert one.full_name == "Test organism STRAIN 1"
    assert one.scaffolds["scaffold_1"].hits == [hits[0], hits[2]]
    assert one.scaffolds["scaffold_2"].hits == [hits[1]]

    assert two.full_name == "Test organism STRAIN 2"
    assert two.scaffolds["scaffold_1"].hits == [hits[3]]
