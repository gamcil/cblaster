"""Tests for genome_parsers.py"""

import pytest

from pathlib import Path

from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from gffutils import Feature

from cblaster import genome_parsers as gp


TEST_DIR = Path(__file__).resolve().parent


@pytest.mark.parametrize(
    "fstart, fend, locations, result",
    [
        (150, 170, [(1, 100), (101, 200), (201, 300)], 1),
        (301, 310, [(1, 100), (101, 200), (201, 300)], None)
    ]
)
def test_find_overlapping_location(fstart, fend, locations, result):
    fl = FeatureLocation(start=fstart, end=fend, strand=1)
    sf = SeqFeature(fl)
    assert gp.find_overlapping_location(sf, locations) == result


@pytest.mark.parametrize(
    "qualifiers, result",
    [
        (dict(), "N.A."),
        (dict(locus_tag="Locus tag"), "Locus tag"),
        (123, TypeError)
    ]
)
def test_find_gene_name(qualifiers, result):
    if type(result) == type and issubclass(result, Exception):
        with pytest.raises(result):
            gp.find_gene_name(qualifiers)
    else:
        assert gp.find_gene_name(qualifiers) == result


def test_find_fasta(tmpdir):
    gff = tmpdir / "test.gff"
    assert gp.find_fasta(gff) is None

    bad = tmpdir / "test.blah"
    bad.write("test")
    assert gp.find_fasta(gff) is None

    fasta = tmpdir / "test.fasta"
    fasta.write("test")
    assert gp.find_fasta(gff) == fasta


def test_parse_fasta_str():
    # Invalid string
    assert gp.parse_fasta_str("abc") == []

    # Actual test case
    fasta = ">Testing\nACGTACGT"
    srs = gp.parse_fasta_str(fasta)
    exp = SeqRecord("ACGTACGT", id="Testing")
    assert exp.id == srs[0].id
    assert exp.seq == srs[0].seq


def test_parse_fasta(tmpdir):
    exp = SeqRecord("ACGTACGT", id="Testing")

    # Does not exist
    with pytest.raises(FileNotFoundError):
        gp.parse_fasta(tmpdir / "doesntexist")

    # Invalid FASTA
    bad = tmpdir / "fasta.fa"
    bad.write("invalid")
    assert gp.parse_fasta(bad) == []

    # Valid FASTA
    good = tmpdir / "fasta.fa"
    good.write(">Testing\nACGTACGT")
    srs = gp.parse_fasta(good)
    assert exp.id == srs[0].id
    assert exp.seq == srs[0].seq


@pytest.mark.parametrize(
    "directives, result",
    [
        ([], {}),
        (["sequence-region a 1 2"], {"a": (1, 2)}),
    ]
)
def test_find_regions(directives, result):
    """Looks for ##sequence-region directives in a list of GFF3 directives."""
    assert gp.find_regions(directives) == result


def test_find_files(tmpdir):
    """find_files correctly parses directory for valid file paths"""
    # Set up directory structure
    # Folder with 2 valid files, 1 invalid
    one = tmpdir / "1"
    one.mkdir()
    gbk = one / "file1.gbk"
    gbk.write("test")
    gff = one / "file2.gff"
    gff.write("test")
    (one / "file2.invalid").write("test")

    # Empty folder
    two = tmpdir / "2"
    two.mkdir()

    # Valid file in root of given path
    base = tmpdir / "file4.gbk"
    base.write("test")

    assert gp.find_files([tmpdir]) == [gbk, gff, base]
    assert gp.find_files([tmpdir], recurse=False) == [base]


def test_parse_file_valid(tmpdir):
    """Files with valid extensions parsed correctly"""
    paths = [
        tmpdir / "test.gbk",
        tmpdir / "test.fasta",
        tmpdir / "test.embl"
    ]
    for path in paths:
        path.write("content")
        assert gp.parse_file(path) == dict(name="test", records=[])


def test_parse_file_gff(tmpdir):
    """GFF files are handled correctly with/without FASTA files"""
    gff = tmpdir / "test.gff"
    gff.write("content")
    with pytest.raises(FileNotFoundError):
        gp.parse_file(gff)

    (tmpdir / "test.fasta").write("content")
    assert gp.parse_file(gff) == dict(name="test", records=[])


def test_parse_file_invalid_ext(tmpdir):
    """parse_file throws ValueError on invalid extensions"""
    invalid = tmpdir / "test.invalid"
    invalid.write("content")
    with pytest.raises(ValueError):
        gp.parse_file(invalid)


def test_parse_cds_features():
    """GFFUtils Feature objects are converted to SeqFeatures and normalised by record
    start"""

    features = [
        Feature(featuretype="gene", start=100, end=1000),
        Feature(featuretype="CDS", start=210, end=800),
    ]
    cds, gene = gp.parse_cds_features(features, 49)

    exp_gene = SeqFeature(FeatureLocation(50, 951, strand=0), type="gene")
    assert gene[0].type == exp_gene.type
    assert gene[0].location.start == exp_gene.location.start
    assert gene[0].location.end == exp_gene.location.end

    exp_cds = SeqFeature(FeatureLocation(160, 751, strand=0), type="CDS")
    assert cds[0].type == exp_cds.type
    assert cds[0].location.start == exp_cds.location.start
    assert cds[0].location.end == exp_cds.location.end


def test_merge_cds_features():
    # Define some feature locations
    exon_one = FeatureLocation(20, 150, strand=1)
    exon_two = FeatureLocation(200, 600, strand=1)
    exon_thr = FeatureLocation(650, 980, strand=1)
    exon_fou = FeatureLocation(1200, 1500, strand=-1)
    exon_fiv = FeatureLocation(1600, 1800, strand=-1)

    # Need ID qualifiers for checking same features
    quals_one = dict(ID=["one"])
    quals_two = dict(ID=["two"])

    # Other SeqFeature objects should merge into these two
    merged_one = SeqFeature(exon_one, type="CDS", qualifiers=quals_one)
    merged_two = SeqFeature(exon_fou, type="CDS", qualifiers=quals_two)

    # Feature list contains other SeqFeatures that will merge into one/two
    features = [
        merged_one,
        SeqFeature(exon_two, type="CDS", qualifiers=quals_one),
        SeqFeature(exon_thr, type="CDS", qualifiers=quals_one),
        merged_two,
        SeqFeature(exon_fiv, type="CDS", qualifiers=quals_two)
    ]

    # Should only have one/two in result, exons should match
    expected = [merged_one, merged_two]
    result = gp.merge_cds_features(features)
    assert result == expected
    assert result[0].location == exon_one + exon_two + exon_thr
    assert result[1].location == exon_fiv + exon_fou  # reversed
