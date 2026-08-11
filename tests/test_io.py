"""Tests for I/O utilities."""

import pytest
import pandas as pd

from optitype.io.data import get_data_path, load_reference_data, get_reference_fasta
from optitype.io import readers


def test_get_data_path():
    """Test that data path can be found."""
    data_path = get_data_path()
    assert data_path.exists()
    assert data_path.is_dir()


def test_load_reference_data():
    """Test loading reference data from Parquet files."""
    table, features = load_reference_data()

    assert isinstance(table, pd.DataFrame)
    assert isinstance(features, pd.DataFrame)

    # Check table structure
    assert "id" in table.columns
    assert "4digit" in table.columns
    assert "locus" in table.columns
    assert len(table) > 0

    # Check features structure
    assert "id" in features.columns
    assert "feature" in features.columns
    assert "number" in features.columns
    assert len(features) > 0


def test_get_reference_fasta():
    """Test getting reference FASTA paths."""
    dna_ref = get_reference_fasta("dna")
    assert dna_ref.exists()
    assert "dna" in dna_ref.name

    rna_ref = get_reference_fasta("rna")
    assert rna_ref.exists()
    assert "rna" in rna_ref.name


def test_get_reference_fasta_invalid():
    """Test that invalid seq_type raises ValueError."""
    with pytest.raises(ValueError):
        get_reference_fasta("invalid")


def test_pysam_to_dataframe_empty_alignments(monkeypatch):
    """Empty BAM/SAM input should return empty DataFrames with valid columns."""

    class DummySam:
        header = {"PG": [{"ID": "yara", "CL": ""}]}
        nreferences = 3
        references = ["HLA:A*01:01", "HLA:B*07:02", "HLA:C*07:02"]

        def __iter__(self):
            return iter(())

    class DummyPysam:
        @staticmethod
        def AlignmentFile(_samfile, _mode):
            return DummySam()

    monkeypatch.setattr(readers, "PYSAM_AVAILABLE", True)
    monkeypatch.setattr(readers, "pysam", DummyPysam(), raising=False)

    pos_df, details_df = readers.pysam_to_dataframe("empty.bam")

    assert pos_df.empty
    assert list(pos_df.columns) == ["HLA:A*01:01", "HLA:B*07:02", "HLA:C*07:02"]
    assert details_df.empty
    assert list(details_df.columns) == ["mismatches", "read_length"]
