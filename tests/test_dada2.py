import pandas as pd
import pytest

from asv_explorer.io.dada2 import (
    load_dada2_counts,
    load_dada2_taxonomy,
    find_zero_count_samples,
    validate_asv_alignment,
)

def test_valid_counts(tmp_path):

    counts = pd.DataFrame({
        "ASV_ID": ["ASV_1", "ASV_2"],
        "Sample_A": [10, 20],
        "Sample_B": [5, 0],
    })

    input_file = tmp_path / "counts.tsv"
    counts.to_csv(input_file, sep="\t", index=False)

    result = load_dada2_counts(input_file)

    assert result.shape == (2, 3)
    assert result["ASV_ID"].is_unique

def test_negative_counts(tmp_path):

    counts = pd.DataFrame({
        "ASV_ID": ["ASV_1", "ASV_2"],
        "Sample_A": [10, -5],
    })

    input_file = tmp_path / "negative.tsv"
    counts.to_csv(input_file, sep="\t", index=False)

    with pytest.raises(ValueError, match="Negative ASV counts"):
        load_dada2_counts(input_file)

def test_fractional_counts(tmp_path):

    counts = pd.DataFrame({
        "ASV_ID": ["ASV_1", "ASV_2"],
        "Sample_A": [10, 2.5],
    })

    input_file = tmp_path / "fractional.tsv"
    counts.to_csv(input_file, sep="\t", index=False)

    with pytest.raises(ValueError, match="Fractional ASV counts detected"):
        load_dada2_counts(input_file)

@pytest.mark.parametrize(
    ("value", "expected_error"),
    [
        ("abc", "Non-numeric ASV counts"),
        (float("inf"), "Non-finite ASV counts"),
        (float("-inf"), "Non-finite ASV counts"),
    ],
)
def test_invalid_numeric_values(tmp_path, value, expected_error):

    counts = pd.DataFrame({
        "ASV_ID": ["ASV_1", "ASV_2"],
        "Sample_A": [10, value],
    })

    input_file = tmp_path / "invalid.tsv"
    counts.to_csv(input_file, sep="\t", index=False)

    with pytest.raises(ValueError, match=expected_error):
        load_dada2_counts(input_file)


def test_zero_count_sample(tmp_path):

    counts = pd.DataFrame({
        "ASV_ID": ["ASV_1", "ASV_2"],
        "Sample_A": [10, 20],
        "Sample_Zero": [0, 0],
    })

    input_file = tmp_path / "zero_sample.tsv"
    counts.to_csv(input_file, sep="\t", index=False)

    result = load_dada2_counts(input_file)

    assert result.shape == (2, 3)
    assert result["Sample_Zero"].sum() == 0

    from asv_explorer.io.dada2 import find_zero_count_samples

    zero_samples = find_zero_count_samples(result)

    assert zero_samples == ["Sample_Zero"]


def test_missing_sample_columns(tmp_path):

    counts = pd.DataFrame({
        "ASV_ID": ["ASV_1", "ASV_2"],
    })

    input_file = tmp_path / "missing_samples.tsv"
    counts.to_csv(input_file, sep="\t", index=False)

    with pytest.raises(ValueError, match="No sample columns detected"):
        load_dada2_counts(input_file)


def test_duplicate_asv_ids(tmp_path):

    counts = pd.DataFrame({
        "ASV_ID": ["ASV_1", "ASV_1"],
        "Sample_A": [10, 20],
    })

    input_file = tmp_path / "duplicates.tsv"
    counts.to_csv(input_file, sep="\t", index=False)

    with pytest.raises(ValueError, match="Duplicated ASV identifiers"):
        load_dada2_counts(input_file)

def test_valid_taxonomy(tmp_path):
    taxonomy = pd.DataFrame({
        "ASV_ID": ["ASV_1", "ASV_2"],
        "Kingdom": ["Bacteria", "Bacteria"],
        "Phylum": ["Proteobacteria", "Bacteroidota"],
    })

    input_file = tmp_path / "taxonomy.tsv"
    taxonomy.to_csv(input_file, sep="\t", index=False)

    result = load_dada2_taxonomy(input_file)

    assert result.shape == (2, 3)
    assert result["ASV_ID"].is_unique


def test_missing_taxonomy_alignment():
    counts = pd.DataFrame({
        "ASV_ID": ["ASV_1", "ASV_2"],
        "Sample_A": [10, 20],
    })

    taxonomy = pd.DataFrame({
        "ASV_ID": ["ASV_1"],
        "Kingdom": ["Bacteria"],
    })

    with pytest.raises(ValueError, match="Missing taxonomy for 1 ASVs"):
        validate_asv_alignment(counts, taxonomy)


def test_duplicate_taxonomy_ids(tmp_path):
    taxonomy = pd.DataFrame({
        "ASV_ID": ["ASV_1", "ASV_1"],
        "Kingdom": ["Bacteria", "Bacteria"],
    })

    input_file = tmp_path / "taxonomy_duplicates.tsv"
    taxonomy.to_csv(input_file, sep="\t", index=False)

    with pytest.raises(
        ValueError,
        match="Duplicated taxonomy ASV identifiers",
    ):
        load_dada2_taxonomy(input_file)