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


def test_find_fasta(tmp_path):
    gff = tmp_path / "test.gff"
    assert gp.find_fasta(gff) is None

    bad = tmp_path / "test.blah"
    bad.write_text("test")
    assert gp.find_fasta(gff) is None

    fasta = tmp_path / "test.fasta"
    fasta.write_text("test")
    assert gp.find_fasta(gff) == fasta


def test_parse_fasta(tmp_path):
    exp = SeqRecord("ACGTACGT", id="Testing")

    # Does not exist
    with pytest.raises(FileNotFoundError):
        gp.parse_fasta(tmp_path / "doesntexist")

    # Invalid FASTA
    bad = tmp_path / "fasta.fa"
    bad.write_text("invalid")
    assert gp.parse_fasta(bad) == []

    # Valid FASTA
    good = tmp_path / "fasta.fa"
    good.write_text(">Testing\nACGTACGT")
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


def test_find_files(tmp_path):
    """find_files correctly parses directory for valid file paths"""
    directory = tmp_path / "folder"
    directory.mkdir()

    # Set up directory structure
    # Folder with 2 valid files, 1 invalid
    one = directory / "1"
    one.mkdir()
    gbk = one / "file1.gbk"
    gbk.write_text("test")
    gff = one / "file2.gff"
    gff.write_text("test")
    (one / "file2.invalid").write_text("test")

    # Empty folder
    two = directory / "2"
    two.mkdir()

    # Valid file in root of given path
    base = directory / "file4.gbk"
    base.write_text("test")

    recursed = gp.find_files([directory], recurse=True)
    assert all(f in recursed for f in [base, gbk, gff])
    assert gp.find_files([directory], recurse=False) == [base]


def test_parse_file_valid(tmp_path):
    """Files with valid extensions parsed correctly"""
    paths = [
        tmp_path / "test.gbk",
        tmp_path / "test.fasta",
        tmp_path / "test.embl"
    ]
    for path in paths:
        path.write_text("content")
        assert gp.parse_file(path) == dict(name="test", records=[])


def test_parse_file_gff(tmp_path):
    """GFF files are handled correctly with/without FASTA files"""
    gff = tmp_path / "test.gff"
    gff.write_text("content")
    with pytest.raises(FileNotFoundError):
        gp.parse_file(gff)

    (tmp_path / "test.fasta").write_text("content")
    assert gp.parse_file(gff) == dict(name="test", records=[])


def test_parse_file_invalid_ext(tmp_path):
    """parse_file throws ValueError on invalid extensions"""
    invalid = tmp_path / "test.invalid"
    invalid.write_text("content")
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
