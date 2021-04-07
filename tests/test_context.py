#!/usr/bin/env python3

"""
Test suite for context.py
"""

import pytest

from pathlib import Path

import io
import requests
import requests_mock

from cblaster import classes, context


TEST_DIR = Path(__file__).resolve().parent

"""
1	INSDC	scaffold_1	14044	14641	-	s1	Prot1	Test organism	STRAIN 1	
2	INSDC	scaffold_2	11815	13459	+	s2	Prot2	Test organism	STRAIN 1	
2	INSDC		11815	13459	+	s2	Prot2	Test organism	STRAIN 123	
3	INSDC	scaffold_1	9656	11184	+	s3	Prot3	Test organism	STRAIN 1	
4	INSDC	scaffold_7	9656	11184	+	s4	Prot4	Test organism STRAIN 2	STRAIN 2	
4	INSDC	scaffold_2	1234	5678	+	s5	Prot4	Test organism	STRAIN 1	
"""


@pytest.fixture()
def hits():
    return [
        classes.Hit("q1", "s1", "1", "1", "0", "1"),
        classes.Hit("q2", "s2", "1", "1", "0", "1"),
        classes.Hit("q3", "s3", "1", "1", "0", "1"),
        classes.Hit("q4", "s4", "1", "1", "0", "1"),
        classes.Hit("q5", "s5", "1", "1", "0", "1"),
        classes.Hit("q6", "s4", "1", "1", "0", "1"),
    ]


@pytest.fixture()
def hit_dict(hits):
    return context.group_hits(hits)


@pytest.fixture()
def subjects(hits):
    return [
        classes.Subject(name="s1", hits=[hits[0]], ipg="1", start=0, end=1000, strand="-"),
        classes.Subject(name="s2", hits=[hits[1]], ipg="2", start=2000, end=3000, strand="+"),
        classes.Subject(name="s3", hits=[hits[2]], ipg="3", start=5000, end=6000, strand="+"),
        classes.Subject(name="s4", hits=[hits[3], hits[4], hits[5]], ipg="4", start=9000, end=10000, strand="+"),
        classes.Subject(name="s5", hits=[hits[3], hits[4], hits[5]], ipg="4", start=14000, end=15000, strand="+"),
    ]


@pytest.fixture()
def ipg_table():
    with (TEST_DIR / "ipg_results.txt").open() as fp:
        table = fp.read()
    return table.split("\n")


@pytest.fixture()
def groups(ipg_table):
    return context.parse_IP_groups(ipg_table)


@pytest.mark.parametrize(
    "unique, gap, minimum, percentage, results",
    [
        (0, 999, 1, 0, [[0], [1], [2], [3], [4]]),
        (0, 1000, 1, 0, [[0, 1], [2], [3], [4]]),
        (2, 2000, 3, 0, [[0, 1, 2]]),
        (4, 2000, 0, 0, []),
        (6, 2000, 0, 0, []),
        (0, 1000, 0, 100, [[0, 1]]),
    ],
)
def test_find_clusters(subjects, unique, gap, minimum, percentage, results):
    groups = context.find_clusters(
        subjects,
        queries=["q1", "q2"],
        unique=unique,
        min_hits=minimum,
        gap=gap,
        percentage=percentage
    )
    for group, result in zip(groups, results):
        assert group == [subjects[i] for i in result]


@pytest.mark.parametrize("unique, gap", [(-1, 100), (1, -1)])
def test_find_clusters_negative_input(subjects, unique, gap):
    with pytest.raises(ValueError):
        x = list(context.find_clusters(subjects, unique=unique, gap=gap))


def test_efetch_IPGs_no_ids():
    with pytest.raises(ValueError):
        context.efetch_IPGs([])


def test_efetch_IPGs_IOError(hits, monkeypatch):
    def mock_ioerror(*args, **kwargs):
        raise IOError
    from Bio import Entrez
    monkeypatch.setattr(Entrez, "efetch", mock_ioerror)
    with pytest.raises(IOError):
        context.efetch_IPGs([hit.subject for hit in hits])


def test_efetch_IPGs_chunks(monkeypatch):
    def mock_chunks(
        _,
        rettype=None,
        retmode=None,
        retmax=None,
        id=None,
    ):
        text = "\n".join(id).encode("utf-8")
        return io.BytesIO(text)
    from Bio import Entrez
    monkeypatch.setattr(Entrez, "efetch", mock_chunks)
    for num in [500, 10001, 1]:
        ids = ["id"] * num
        table = context.efetch_IPGs(ids)
        assert len(table) == num


def test_efetch_IPGs_output(tmp_path, monkeypatch):
    def mock_output(*args, **kwargs):
        text = "testing".encode("utf-8")
        return io.BytesIO(text)
    from Bio import Entrez
    monkeypatch.setattr(Entrez, "efetch", mock_output)
    test_out = tmp_path / "out.tsv"
    context.efetch_IPGs(["placeholder"], output_file=test_out)
    assert test_out.read_text() == "testing"


def test_parse_IPG_table(hits, subjects):
    results = TEST_DIR / "ipg_results.txt"

    with open(results) as handle:
        parsed = context.parse_IPG_table(handle, hits)

    assert len(parsed) == 2, "Expected 2 organism objects"
    one, two = parsed

    assert one.name == two.name == "Test organism"
    assert one.strain == "STRAIN 1"
    assert two.strain == "STRAIN 2"

    assert len(one.scaffolds) == 2, "Expected 2 scaffolds in organism 1"
    assert len(two.scaffolds) == 1, "Expected 1 scaffold in organism 2"

    assert one.full_name == "Test organism STRAIN 1"

    assert len(one.scaffolds["scaffold_1"].subjects) == 2, "Expects 2 subjects on scaf 1"
    assert one.scaffolds["scaffold_1"].subjects[0].name == subjects[0].name
    assert one.scaffolds["scaffold_1"].subjects[1].name == subjects[2].name

    assert len(one.scaffolds["scaffold_2"].subjects) == 2, "Expects 2 subjects on scaf 2"
    assert one.scaffolds["scaffold_2"].subjects[0].name == subjects[1].name
    assert one.scaffolds["scaffold_2"].subjects[1].name == subjects[4].name

    assert two.full_name == "Test organism STRAIN 2"
    assert two.scaffolds["scaffold_7"].subjects[0].name == subjects[3].name


def test_parse_IP_groups(ipg_table):
    x = context.parse_IP_groups(ipg_table)
    assert len(x) == 4, "Expected 4 groups"
    assert len(x["2"]) == 2, "Group 2 has 2 rows"


@pytest.mark.parametrize(
    "group, length",
    [("1", 1), ("2", 1), ("3", 1), ("4", 3)],
)
def test_find_IPG_hits(groups, hits, hit_dict, group, length):
    x = context.find_IPG_hits(groups[group], hit_dict)
    assert len(x) == length, "Hit group length mismatch"
