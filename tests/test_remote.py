#!/usr/bin/env python3

"""
Test suite for remote module
"""

import pytest
import requests_mock

from pathlib import Path


from cblaster import remote, helpers


TEST_DIR = Path(__file__).resolve().parent


def test_start_no_input():
    with pytest.raises(ValueError):
        # No query_file/ids
        remote.start()


@pytest.fixture()
def start_response():
    return (TEST_DIR / "start_response.html").read_text()


def test_start(start_response, monkeypatch):
    def mock_sequences(query_file, query_ids):
        return {'seq1': 'TEST', 'seq2': 'TEST'}

    monkeypatch.setattr(helpers, "get_sequences", mock_sequences)

    with requests_mock.Mocker() as mock:
        mock.post(remote.BLAST_API_URL, text=start_response)

        # Ensure RID/RTOE is returned
        result = remote.start(
            query_ids=["seq1", "seq2"],
            entrez_query="Aspergillus[ORGN]"
        )

        assert result == ("VCZM3MWB014", 18)

        # Check correct request URL
        assert mock.request_history[0].url == (
            "https://blast.ncbi.nlm.nih.gov/Blast.cgi?"
            "CMD=PUT"
            "&DATABASE=nr"
            "&PROGRAM=blastp"
            "&FILTER=F"
            "&EXPECT=0.1"
            "&GAPCOSTS=11+1"
            "&MATRIX=BLOSUM62"
            "&HITLIST_SIZE=500"
            "&ALIGNMENTS=500"
            "&DESCRIPTIONS=500"
            "&WORD_SIZE=6"
            "&COMPOSITION_BASED_STATISTICS=2"
            "&ENTREZ_QUERY=Aspergillus%5BORGN%5D"
            "&THRESHOLD=11"
        )


def test_start_blastn_options(start_response, monkeypatch):
    def mock_sequences(query_file, query_ids):
        return {'seq1': 'TEST', 'seq2': 'TEST'}

    monkeypatch.setattr(helpers, "get_sequences", mock_sequences)

    with requests_mock.Mocker() as mock:
        mock.post(remote.BLAST_API_URL, text=start_response)

        # megablast, nucl_* are blastn options, threshold is only BLASTp
        remote.start(
            query_ids=["seq1"],
            program="blastn",
            megablast=True,
            nucl_penalty=99,
            nucl_reward=99,
            threshold=99,
        )

        # Check correct request URL
        request = mock.request_history[0]
        assert "THRESHOLD" not in request.url  # Only blastp
        assert all(
            part in request.url
            for part in ["NUCL_PENALTY=99", "NUCL_REWARD=99", "MEGABLAST=on"]
        )


@pytest.fixture()
def check_response():
    return (TEST_DIR / "check_response.html").read_text()


def test_check(check_response):
    with requests_mock.Mocker() as mock:
        mock.get(remote.BLAST_API_URL, text=check_response)

        # Finds Status=READY and ThereAreHits=yes
        assert remote.check("VCZM3MWB014") is True

        # Check correct request URL
        assert mock.request_history[0].url == (
            "https://blast.ncbi.nlm.nih.gov/Blast.cgi?"
            "CMD=Get"
            "&RID=VCZM3MWB014"
            "&FORMAT_OBJECT=SearchInfo"
        )


@pytest.mark.parametrize(
    "text", ["Status=UNKNOWN\n", "Status=FAILED\n", "Status=READY\nThereAreHits=no\n"]
)
def test_check_failed(text):
    with requests_mock.Mocker() as mock, pytest.raises(ValueError):
        mock.get(remote.BLAST_API_URL, text=text)
        remote.check("RID")


def test_check_waiting():
    with requests_mock.Mocker() as mock:
        mock.get(remote.BLAST_API_URL, text="Status=WAITING\n")
        assert remote.check("RID") is False


@pytest.fixture()
def retrieve_response():
    return (TEST_DIR / "retrieve_response.html").read_text()


def test_retrieve(retrieve_response):
    with requests_mock.Mocker() as mock:
        mock.get(remote.BLAST_API_URL, text=retrieve_response)

        result = remote.retrieve("RID")

        # Make sure we've removed non-TSV cruft
        assert len(result) == 300
        assert not any(row.startswith(("#", "<", " ", "Qblast", "-")) for row in result)
        assert mock.request_history[0].url == (
            "https://blast.ncbi.nlm.nih.gov/Blast.cgi?"
            "CMD=Get"
            "&RID=RID"
            "&FORMAT_TYPE=Tabular"
            "&FORMAT_OBJECT=Alignment"
            "&HITLIST_SIZE=5000"
            "&ALIGNMENTS=5000"
            "&DESCRIPTIONS=5000"
            "&NCBI_GI=F"
        )


def test_poll_success(monkeypatch):
    def patch_check(rid):
        return True

    monkeypatch.setattr(remote, "check", patch_check)

    assert remote.check("RID") is True


def test_poll_retry_limit(monkeypatch):
    def returns_false(rid):
        return False

    monkeypatch.setattr(remote, "check", returns_false)

    with pytest.raises(ValueError):
        remote.poll("RID", delay=0, max_retries=2)


@pytest.fixture
def query_file():
    return TEST_DIR / "test.faa"


def test_parse_empty_handle(query_file):
    with pytest.raises(ValueError):
        remote.parse([], query_file=query_file)


def test_parse(query_file):
    # length of QBE85648 == 179
    result = [
        # qid sid pid len mismatch gapopen qstart qend sstart ssend evalue bitscore
        "QBE85648.1\tHIT1\t100.000\t179\t0\t0\t1\t179\t1\t179\t1.38e-127\t365\t100.00",
        "QBE85648.1\tHIT2\t20.000\t179\t0\t0\t1\t179\t1\t179\t1.38e-127\t365\t100.00",
        "QBE85648.1\tHIT3\t100.000\t179\t0\t0\t150\t179\t1\t179\t1.38e-127\t365\t100.00",
        "QBE85648.1\tHIT4\t100.000\t179\t0\t0\t1\t179\t1\t179\t0.011\t365\t100.00",
    ]

    hits = remote.parse(result, query_file=query_file)

    # Default thresholds are 30% identity, 50% coverage, 0.01 evalue
    # so only the first hit should be saved
    assert len(hits) == 1
    assert hits[0].query == "QBE85648.1"
    assert hits[0].subject == "HIT1"
    assert hits[0].identity == 100.0
    assert hits[0].coverage == 100.0
    assert hits[0].bitscore == 365.0
    assert hits[0].evalue == 1.38e-127
