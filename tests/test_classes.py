#!/usr/bin/env python3

"""
Test suite for classes.
"""

import pytest

from cblaster import classes


@pytest.fixture
def hits():
    hits = [
        classes.Hit("q1", "s1", "70.90", "54.60", "0.0", "500.30", 100, 1000, "+"),
        classes.Hit("q1", "s2", "70.90", "54.60", "1.06e-29", "100.30", 100, 1000, "+"),
        classes.Hit("q2", "s3", "70.90", "54.60", "1.02e-30", "200.30", 100, 1000, "+"),
        classes.Hit("q3", "s5", "70.90", "54.60", "0.01", "500.30", 100, 1000, "+"),
        classes.Hit("q4", "s1", "70.90", "54.60", "0.005", "500.30", 100, 1000, "+"),
    ]

    return hits


@pytest.mark.parametrize(
    "text, symbol, result",
    [("test", "-", "test\n----"), ("test2", "=", "test2\n=====")],
)
def test_generate_header_string(text, symbol, result):
    assert classes.generate_header_string(text, symbol) == result


def test_generate_cluster_table(hits):
    assert classes.generate_cluster_table(hits, decimals=4, show_headers=True) == (
        "Query  Subject  Identity  Coverage  E-value   Bitscore  Start  End   Strand\n"
        "q1     s1       70.9      54.6      0         500.3     100    1000  +     \n"
        "q1     s2       70.9      54.6      1.06e-29  100.3     100    1000  +     \n"
        "q2     s3       70.9      54.6      1.02e-30  200.3     100    1000  +     \n"
        "q3     s5       70.9      54.6      0.01      500.3     100    1000  +     \n"
        "q4     s1       70.9      54.6      0.005     500.3     100    1000  +     "
    )


@pytest.mark.parametrize(
    "index, decimals, result",
    [
        (0, 4, ["q1", "s1", "70.9", "54.6", "0", "500.3", "100", "1000", "+"]),
        (0, 0, ["q1", "s1", "71", "55", "0", "500", "100", "1000", "+"]),
        (1, 2, ["q1", "s2", "70.9", "54.6", "1.1e-29", "100.3", "100", "1000", "+"]),
    ],
)
def test_hit_values(hits, index, decimals, result):
    assert hits[index].values(decimals) == result


def test_hit_str(hits):
    assert str(hits[0]) == "HIT: q1 - s1 [70.9, 54.6]"


def test_hit_instantiation(hits):
    assert [str, str, float, float, float, float, int, int, str] == [
        type(getattr(hits[0], val))
        for val in [
            "query",
            "subject",
            "identity",
            "coverage",
            "evalue",
            "bitscore",
            "start",
            "end",
            "strand",
        ]
    ]


@pytest.fixture()
def scaffold(hits):
    return classes.Scaffold("test_scaffold", hits)


def test_scaffold_instantation(scaffold, hits):
    assert scaffold.accession == "test_scaffold"
    assert scaffold.hits == hits
    assert scaffold.clusters == []


def test_scaffold_str(scaffold):
    assert str(scaffold) == "SCAFFOLD: test_scaffold [5 hits in 0 clusters]"


SCAFFOLD_SUMMARY = (
    "test_scaffold\n"
    "-------------\n"
    "Query  Subject  Identity  Coverage  E-value   Bitscore  Start  End   Strand\n"
    "q1     s1       70.9      54.6      0         500.3     100    1000  +     \n"
    "q1     s2       70.9      54.6      1.06e-29  100.3     100    1000  +     \n"
    "q2     s3       70.9      54.6      1.02e-30  200.3     100    1000  +     \n"
    "\n"
    "Query  Subject  Identity  Coverage  E-value  Bitscore  Start  End   Strand\n"
    "q3     s5       70.9      54.6      0.01     500.3     100    1000  +     \n"
    "q4     s1       70.9      54.6      0.005    500.3     100    1000  +     "
)


def test_scaffold_summary(scaffold):
    with pytest.raises(ValueError):
        # No clusters
        scaffold.summary()

    scaffold.clusters.append(scaffold.hits[:3])
    scaffold.clusters.append(scaffold.hits[3:])

    assert scaffold.summary() == SCAFFOLD_SUMMARY


@pytest.fixture()
def organism():
    return classes.Organism("test_organism", "TEST 123")


def test_organism_full_name(organism):
    assert organism.full_name == "test_organism TEST 123"

    # Test again if strain already in organism name, as some are in NCBI
    organism.name = "test_organism TEST 123"
    assert organism.full_name == "test_organism TEST 123"


def test_organism_instantation(organism):
    assert organism.name == "test_organism"
    assert organism.strain == "TEST 123"
    assert organism.scaffolds == {}


def test_organism_str_no_scaffolds(organism):
    assert str(organism) == "ORGANISM: test_organism TEST 123 [0 hits on 0 scaffolds]"


@pytest.fixture()
def org_with_clusters(scaffold):
    organism = classes.Organism("test_organism", "TEST 123")
    organism.scaffolds[scaffold.accession] = scaffold
    organism.scaffolds["test_scaffold"].clusters = [
        organism.scaffolds["test_scaffold"].hits[:3],
        organism.scaffolds["test_scaffold"].hits[3:],
    ]
    return organism


def test_organism_str_with_scaffolds(org_with_clusters):
    assert (
        str(org_with_clusters)
        == "ORGANISM: test_organism TEST 123 [5 hits on 1 scaffolds]"
    )


def test_organism_count_hit_clusters(organism, org_with_clusters):
    assert organism.count_hit_clusters() == 0
    assert org_with_clusters.count_hit_clusters() == 2


def test_organism_summary(organism, org_with_clusters):
    with pytest.raises(ValueError):
        organism.summary()

    assert org_with_clusters.summary() == (
        "test_organism TEST 123\n======================\n" + SCAFFOLD_SUMMARY
    )
