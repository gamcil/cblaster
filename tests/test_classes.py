# !/usr/bin/env python3

"""
Test suite for classes.
"""

import pytest
import difflib

from abc import ABC
from pathlib import Path
from tempfile import TemporaryFile
from cblaster import classes


TEST_DIR = Path(__file__).resolve().parent
DATA_DIR = TEST_DIR / "data"


# HIT TESTS
@pytest.fixture
def hits():
    hits = [
        classes.Hit("q1", "s1", "70.90", "54.60", "0.0", "500.30"),
        classes.Hit("q1", "s2", "5.90", "21.60", "1.06e-29", "100.30"),
        classes.Hit("q2", "s3", "2.90", "34.60", "1.02e-30", "200.30"),
        classes.Hit("q3", "s5", "1123.90", "12.60", "0.01", "500.30"),
        classes.Hit("q4", "s1", "0.90", "1.60", "0.005", "500.30"),
        classes.Hit("q1", "s1", "70.90", "54.60", "0.0", "500.30"),
    ]

    return hits


@pytest.fixture
def subjects(hits):
    return [
        classes.Subject("id1", hits[:2], "subject1", "ipg1", "0", "100", -1),
        classes.Subject("id1", [hits[1], hits[0]], "subject1", "ipg1", 0, 100, -1),
        classes.Subject("id2", hits[2:3], "subject2", "ipg2", 1200, 6000, 1),
        classes.Subject("id3", hits[3:4], "subject3", "ipg3", 500, "1000", 1),
        classes.Subject("id4", hits[4:5], "subject4", "ipg4", 6453, 8000, 1),
        classes.Subject("id5", name="intermediate subject", start=8100, end=9000),
    ]


@pytest.fixture
def clusters(subjects):
    return [
        classes.Cluster([0, 1], subjects[:2], number=1),
        classes.Cluster([0, 1], subjects[:2], score=100, start=50, end=200, number=2),
        classes.Cluster([0, 2, 3, 4], [subjects[0], *subjects[2:5]], number=3),
        classes.Cluster([0, 1], [subjects[1], subjects[0]], number=4),
        classes.Cluster(
            [0, 1],
            [subjects[1], subjects[0]],
            intermediate_genes=[subjects[5]],
            number=5,
        ),
    ]


@pytest.fixture
def scaffolds(subjects, clusters):
    return [
        classes.Scaffold("scaffold1", [clusters[2], clusters[0]], subjects[:5]),
    ]


@pytest.fixture
def scaffold_summaries():
    with (DATA_DIR / "scaffold_summaries.txt").open() as fp:
        summary_text = fp.read()
    return [s.strip("\n") for s in summary_text.split("\n>>\n") if len(s) > 0]


@pytest.fixture
def organisms(scaffolds):
    return [
        classes.Organism(
            "organism1",
            "strain_organism1",
            {scaffold.accession: scaffold for scaffold in scaffolds[0:1]},
        ),
        classes.Organism(
            None, None, {scaffold.accession: scaffold for scaffold in scaffolds[0:1]}
        ),
        classes.Organism(
            "organism2",
            None,
            {scaffold.accession: scaffold for scaffold in scaffolds[1:2]},
        ),
    ]


@pytest.fixture
def organism_summaries():
    with (DATA_DIR / "organism_summaries.txt").open() as fp:
        summary_text = fp.read()
    return [s.strip("\n") for s in summary_text.split("\n>>\n") if len(s) > 0]


@pytest.fixture
def sessions(organisms, clusters):
    return [
        classes.Session(
            query=clusters[0],
            queries=["q1", "q2"],
            organisms=organisms[0:1],
            # sequences={"q1": "PRQTEINQNE", "q2": "PRQTEINTWQ"},
        ),
        classes.Session(
            query=clusters[1],
            queries=["q1", "q2"],
            organisms=organisms[2:3],
            # sequences={"q1": "SQMESEQVENCE", "q2": "SQMEOTHERSEQVENCE"},
            params={"dummy1": 1},
        ),
    ]


@pytest.fixture
def session_summaries():
    with (DATA_DIR / "session_summaries.txt").open() as fp:
        summary_text = fp.read()
    return [s for s in summary_text.split("\n>>\n") if len(s) > 0]


# TEST HITS
def test_hit_instantiation(hits):
    assert isinstance(hits[0], classes.Serializer)
    assert [str, str, float, float, float, float] == [
        type(getattr(hits[0], val))
        for val in ["query", "subject", "identity", "coverage", "evalue", "bitscore"]
    ]


@pytest.mark.parametrize(
    "index, decimals, result",
    [
        (0, 4, ["q1", "s1", "70.9", "54.6", "0", "500.3"]),
        (0, 0, ["q1", "s1", "71", "55", "0", "500"]),
        (1, 2, ["q1", "s2", "5.9", "21.6", "1.1e-29", "100.3"]),
    ],
)
def test_hit_values(hits, index, decimals, result):
    assert hits[index].values(decimals) == result


def test_hit_str(hits):
    assert str(hits[0]) == "Hit: q1 - s1: 70.9/54.6"


def test_hit_equals(hits):
    with pytest.raises(NotImplementedError):
        hits[0] == "invalid value"
    assert hits[0] == hits[5]
    assert hits[0] != hits[4]


def test_hit_copy(hits):
    copy_hit = hits[0].copy(test="test_value")
    assert ["q1", "s1", 70.9, 54.6, 0.0, 500.3, "test_value"] == [
        getattr(copy_hit, val)
        for val in [
            "query",
            "subject",
            "identity",
            "coverage",
            "evalue",
            "bitscore",
            "test",
        ]
    ]


def test_hit_to_dict(hits):
    assert hits[0].to_dict() == {
        "query": "q1",
        "subject": "s1",
        "identity": 70.9,
        "coverage": 54.6,
        "evalue": 0.0,
        "bitscore": 500.3,
    }


def test_hit_from_dict():
    from_dict_class = classes.Hit.from_dict(
        {
            "query": "q1",
            "subject": "s1",
            "identity": 70.9,
            "coverage": 54.6,
            "evalue": 0.0,
            "bitscore": 500.3,
        }
    )
    assert ["q1", "s1", 70.9, 54.6, 0.0, 500.3] == [
        getattr(from_dict_class, val)
        for val in ["query", "subject", "identity", "coverage", "evalue", "bitscore",]
    ]


# SUBJECT TESTS
def test_subject_instantiation(subjects):
    assert isinstance(subjects[0], classes.Serializer)
    assert [list, str, str, int, int, int] == [
        type(getattr(subjects[0], val))
        for val in ["hits", "ipg", "name", "start", "end", "strand"]
    ]


def test_subject_empty_instantiation():
    assert [[], None, None, None, None, None] == [
        getattr(classes.Subject(), val)
        for val in ["hits", "ipg", "name", "start", "end", "strand"]
    ]


def test_subject_equals(subjects):
    with pytest.raises(NotImplementedError):
        subjects[0] == "invalid value"
    assert subjects[0] == subjects[1]
    assert subjects[0] != subjects[2]


def test_subject_to_dict(hits, subjects):
    assert subjects[3].to_dict() == {
        "id": "id3",
        "hits": [hits[3].to_dict()],
        "name": "subject3",
        "ipg": "ipg3",
        "start": 500,
        "end": 1000,
        "strand": 1,
        "sequence": None,
    }


@pytest.mark.parametrize(
    "index, decimals, result",
    [
        (
            0,
            4,
            [
                ("q1", "s1", "70.9", "54.6", "0", "500.3", "0", "100", "-"),
                ("q1", "s2", "5.9", "21.6", "1.06e-29", "100.3", "0", "100", "-"),
            ],
        ),
        (
            0,
            0,
            [
                ("q1", "s1", "71", "55", "0", "500", "0", "100", "-"),
                ("q1", "s2", "6", "22", "1e-29", "100", "0", "100", "-"),
            ],
        ),
    ],
)
def test_subject_values(subjects, index, decimals, result):
    assert subjects[index].values(decimals) == result


def test_subject_from_dict(hits):
    from_dict_class = classes.Subject.from_dict(
        {
            "hits": [hits[3].to_dict()],
            "name": "subject3",
            "ipg": "ipg3",
            "start": 500,
            "end": 1000,
            "strand": "+",
        }
    )
    assert [[hits[3]], "subject3", "ipg3", 500, 1000, "+"] == [
        getattr(from_dict_class, val)
        for val in ["hits", "name", "ipg", "start", "end", "strand",]
    ]


# CLUSTER TESTS
def test_cluster_instantiation(clusters):
    assert isinstance(clusters[0], classes.Serializer)
    assert [list, list, list, float, int, int, int] == [
        type(getattr(clusters[0], val))
        for val in [
            "indices",
            "subjects",
            "intermediate_genes",
            "score",
            "start",
            "end",
            "number",
        ]
    ]


def test_cluster_empty_instantiation(subjects):
    # other values like subjects are forced to be a number
    assert [[], []] == [
        getattr(classes.Cluster(subjects=[subjects[0]]), val)
        for val in ["indices", "intermediate_genes",]
    ]


def test_cluster_iter(subjects, clusters):
    assert subjects[:2] == [s for s in clusters[0]]


def test_cluster_len(clusters):
    assert len(clusters[0]) == 2
    assert len(clusters[1]) == 2
    assert len(clusters[2]) == 4


def test_cluster_equals(clusters):
    with pytest.raises(NotImplementedError):
        clusters[0] == "invalid value"
    assert clusters[0] != clusters[1]
    assert clusters[0] == clusters[3]


def test_cluster_calculate_score(clusters):
    assert clusters[0].calculate_score() == 2.10006
    assert clusters[2].calculate_score() == 4.17012
    assert clusters[2].calculate_score(["q1", "q3", "q2", "q4"]) == 5.17012
    assert clusters[2].calculate_score(["q1", "q2", "q3", "q4"]) == 7.17012


def test_cluster_to_dict(clusters):
    assert clusters[0].to_dict() == {
        "number": 1,
        "intermediate_genes": [],
        "indices": [0, 1],
        "score": 2.10006,
        "start": 0,
        "end": 100,
        "subjects": []
    }


def test_cluster_from_dict(subjects):
    from_dict_class = classes.Cluster.from_dict(
        {
            "indices": [0, 1],
            "intermediate_genes": [],
            "score": 2.10006,
            "start": 0,
            "end": 100,
            "number": 100,
        },
        subjects=subjects[:2],
    )
    assert [[0, 1], subjects[:2], [], 2.10006, 0, 100, 100] == [
        getattr(from_dict_class, val)
        for val in [
            "indices",
            "subjects",
            "intermediate_genes",
            "score",
            "start",
            "end",
            "number",
        ]
    ]


# TEST SCAFFOLD
def test_scaffold_instantiation(scaffolds):
    assert isinstance(scaffolds[0], classes.Serializer)
    assert [str, list, list] == [
        type(getattr(scaffolds[0], val))
        for val in ["accession", "subjects", "clusters"]
    ]


def test_scaffold_empty_instantiation():
    assert ["test_name", [], []] == [
        getattr(classes.Scaffold("test_name"), val)
        for val in ["accession", "subjects", "clusters"]
    ]


def test_scaffold_string(scaffolds):
    assert str(scaffolds[0]) == "SCAFFOLD: scaffold1 [5 hits in 2 clusters]"


def test_scaffold_add_clusters(subjects, scaffolds):
    scaffolds[0].add_clusters([subjects[1:4]])
    assert len(scaffolds[0].clusters) == 3
    assert all(
        subject
        in subjects[1:4]
        for subject in scaffolds[0].clusters[1].subjects
    )


@pytest.mark.parametrize(
    "index, hide_headers, delimiter, decimals, result_index",
    [
        (0, False, None, 4, 0),
        (0, True, None, 4, 1),
        (0, False, "____", 4, 2),
        (0, False, None, 0, 3),
    ],
)
def test_scaffold_summary(
    scaffolds,
    index,
    hide_headers,
    delimiter,
    decimals,
    result_index,
    scaffold_summaries,
):
    assert (
        scaffolds[index].summary(hide_headers, delimiter, decimals)
        == scaffold_summaries[result_index]
    )


def test_scaffold_to_dict(scaffolds):
    assert scaffolds[0].to_dict() == {
        "accession": "scaffold1",
        "clusters": [cluster.to_dict() for cluster in scaffolds[0].clusters],
        "subjects": [subject.to_dict() for subject in scaffolds[0].subjects],
    }


def test_scaffold_from_dict(subjects, clusters):
    from_dict_class = classes.Scaffold.from_dict(
        {
            "accession": "scaffold2",
            "clusters": [clusters[0].to_dict()],
            "subjects": [subject.to_dict() for subject in subjects[:2]],
        }
    )
    assert ["scaffold2", [clusters[0]], subjects[:2]] == [
        getattr(from_dict_class, val) for val in ["accession", "clusters", "subjects"]
    ]


# TEST ORGANISMS
def test_organism_instantiation(organisms):
    assert isinstance(organisms[0], classes.Serializer)
    assert [str, str, dict] == [
        type(getattr(organisms[0], val)) for val in ["name", "strain", "scaffolds"]
    ]


def test_organism_empty_instantiation():
    assert ["test_name", "test_strain", {}] == [
        getattr(classes.Organism("test_name", "test_strain"), val)
        for val in ["name", "strain", "scaffolds"]
    ]


def test_organism_string(organisms):
    assert (
        str(organisms[0])
        == "ORGANISM: organism1 strain_organism1 [5 subjects on 1 scaffolds]"
    )


def test_organism_clusters(clusters, organisms):
    assert organisms[0].clusters == [clusters[2], clusters[0]]


def test_organism_total_hit_clusters(organisms):
    assert organisms[0].total_hit_clusters == 2


@pytest.mark.parametrize(
    "index, hide_headers, delimiter, decimals, result_index",
    [
        (0, False, None, 4, 0),
        (0, True, None, 4, 1),
        (0, False, "____", 4, 2),
        (0, False, None, 0, 3),
    ],
)
def test_organism_summary(
    organisms,
    index,
    hide_headers,
    delimiter,
    decimals,
    result_index,
    organism_summaries,
):
    assert (
        organisms[index].summary(
            hide_headers=hide_headers, delimiter=delimiter, decimals=decimals
        )
        == organism_summaries[result_index]
    )


def test_organism_full_name(organisms):
    assert organisms[0].full_name == "organism1 strain_organism1"
    assert organisms[1].full_name == "No organism"
    assert organisms[2].full_name == "organism2"


def test_organism_to_dict(organisms):
    assert organisms[0].to_dict() == {
        "name": "organism1",
        "strain": "strain_organism1",
        "scaffolds": [
            scaffold.to_dict() for scaffold in organisms[0].scaffolds.values()
        ],
    }


def test_organism_from_dict(scaffolds):
    from_dict_class = classes.Organism.from_dict(
        {
            "name": "organism1",
            "strain": "strain_organism1",
            "scaffolds": [scaffold.to_dict() for scaffold in scaffolds[0:1]],
        }
    )
    assert ["organism1", "strain_organism1"] == [
        getattr(from_dict_class, val) for val in ["name", "strain",]
    ]
    # make sure that the scaffolds are the same without writing an __eq__ method that would not realy serve a purpose
    from_dict_scaffold = list(from_dict_class.scaffolds.values())[0]
    assert from_dict_scaffold.accession == scaffolds[0].accession
    assert all(
        from_dict_scaffold.subjects[index] == scaffolds[0].subjects[index]
        for index in range(len(from_dict_scaffold.subjects))
    )
    assert all(
        from_dict_scaffold.clusters[index] == scaffolds[0].clusters[index]
        for index in range(len(from_dict_scaffold.clusters))
    )


# TEST SESSIONS
def test_session_instantiation(sessions):
    assert isinstance(sessions[0], classes.Serializer)
    assert [list, dict, list, dict] == [
        type(getattr(sessions[0], val))
        for val in ["queries", "params", "organisms", "sequences",]
    ]


def test_session_empty_instantiation():
    assert [[], {}, [], {}] == [
        getattr(classes.Session(), val)
        for val in ["queries", "params", "organisms", "sequences",]
    ]


def test_session_add(sessions, organisms):
    combined_session_1_2 = sessions[0] + sessions[1]
    assert [["q1", "q2"], {}] == [
        getattr(combined_session_1_2, val)
        for val in ["queries", "params"]
    ]
    # same organism names can assume same organisms in this case
    assert all(
        combined_session_1_2.organisms[index].name == name
        for index, name in enumerate([o.name for o in [organisms[0], organisms[2]]])
    )

    combined_session_2_1 = sessions[1] + sessions[0]
    assert [
        ["q1", "q2"],
        {"dummy1": 1},
    ] == [
        getattr(combined_session_2_1, val)
        for val in ["queries", "params"]
    ]
    # same organism names can assume same organisms in this case
    assert all(
        combined_session_2_1.organisms[index].name == name
        for index, name in enumerate([o.name for o in [organisms[2], organisms[0]]])
    )


def test_session_to_dict(sessions, clusters):
    assert sessions[1].to_dict() == {
        "query": clusters[1].to_dict(save_subjects=True),
        "queries": ["q1", "q2"],
        "params": {"dummy1": 1},
        "organisms": [organism.to_dict() for organism in sessions[1].organisms],
    }


def compare_sessions(one, two):
    """Compares two Session objects"""

    # Compare search attributes
    attributes = ["queries", "params"]
    for attr in attributes:
        assert getattr(one, attr) == getattr(two, attr)

    # Compare organisms
    for a, b in zip(one.organisms, two.organisms):
        compare_organisms(a, b)


def compare_organisms(one, two):
    assert one.full_name == two.full_name, "Name mismatch"
    for scaffold_one, scaffold_two in zip(
        one.scaffolds.values(), two.scaffolds.values()
    ):
        compare_scaffolds(scaffold_one, scaffold_two)


def compare_scaffolds(one, two):
    assert one.accession == two.accession
    for a, b in zip(one.subjects, two.subjects):
        compare_objects(a, b)
    for a, b in zip(one.clusters, two.clusters):
        compare_objects(a, b)


def compare_objects(one, two):
    assert one.to_dict() == two.to_dict()


def test_session_from_file(sessions):
    from_file_session = classes.Session.from_file(DATA_DIR / "test_session1_file.json")
    test_session = sessions[0]
    compare_sessions(from_file_session, test_session)


def test_session_from_files(sessions):
    files = [
        DATA_DIR / "test_session1_file.json",
        DATA_DIR / "test_session2_file.json",
    ]
    from_files_session = classes.Session.from_files(files)
    combined_session = sessions[0] + sessions[1]
    compare_sessions(combined_session, from_files_session)


def test_session_from_dict(sessions, clusters):
    data = {
        "query": clusters[0].to_dict(),
        "queries": ["q1", "q2"],
        "params": {},
        "organisms": [organism.to_dict() for organism in sessions[0].organisms],
    }
    from_dict_session = classes.Session.from_dict(data)
    compare_sessions(sessions[0], from_dict_session)


@pytest.mark.parametrize(
    "index, form, hide_headers, delimiter, decimals, result_index",
    [
        (0, "summary", False, None, 4, 0),
        (0, "summary", True, None, 4, 1),
        (0, "summary", False, "____", 4, 2),
        (0, "summary", False, None, 0, 3),
        (0, "binary", False, None, 4, 4),
        (0, "binary", True, None, 4, 5),
        (0, "binary", False, "____", 4, 6),
        (0, "binary", False, None, 0, 7),
    ],
)
def test_session_format(
    index,
    form,
    hide_headers,
    delimiter,
    decimals,
    result_index,
    sessions,
    session_summaries,
):
    with TemporaryFile(mode="w+") as f:
        sessions[index].format(
            form=form,
            fp=f,
            hide_headers=hide_headers,
            delimiter=delimiter,
            decimals=decimals,
        )
        f.seek(0)
        assert f.read() == session_summaries[result_index]


# TEST SERIALIZER
def test_serializer_is_abstract():
    assert issubclass(classes.Serializer, ABC)
