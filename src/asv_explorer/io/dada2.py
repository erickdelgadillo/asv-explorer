"""Input and validation functions for DADA2 ASV count tables."""

from os import path
from pathlib import Path

from decimal import Decimal, InvalidOperation
import pandas as pd
import csv

def validate_tsv_headers(path: str | Path) -> None:
    """Reject duplicate or empty column names before pandas loads a TSV."""

    path = Path(path)

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        headers = next(reader, [])

    if not headers:
        raise ValueError("Missing TSV header")

    if any(not name.strip() for name in headers):
        raise ValueError("Empty column names detected")

    if len(headers) != len(set(headers)):
        raise ValueError("Duplicated column names detected")

def parse_asv_count(value: str) -> int:
    """Validate and convert an ASV count without floating-point rounding."""

    try:
        number = Decimal(value)
    except InvalidOperation as exc:
        raise ValueError("Non-numeric ASV counts detected") from exc

    if not number.is_finite():
        raise ValueError("Non-finite ASV counts detected")

    if number < 0:
        raise ValueError("Negative ASV counts detected")

    if number != number.to_integral_value():
        raise ValueError("Fractional ASV counts detected")

    max_int64 = 2**63 - 1

    if number > max_int64:
        raise ValueError("ASV counts exceed supported integer range")

    return int(number)

def load_dada2_counts(path: str | Path) -> pd.DataFrame:
    """Load and validate a DADA2 ASV count table."""

    path = Path(path)

    validate_tsv_headers(path)

    counts = pd.read_csv(path, sep="\t", dtype=str)

    # 1. Validate table structure
    if counts.empty:
        raise ValueError("ASV count table is empty")

    if counts.columns.duplicated().any():
        raise ValueError("Duplicated column names detected")

    if "ASV_ID" not in counts.columns:
        raise ValueError("Missing required column: ASV_ID")

    sample_columns = counts.columns.drop("ASV_ID")

    if len(sample_columns) == 0:
        raise ValueError("No sample columns detected")

    if counts["ASV_ID"].isna().any() or counts["ASV_ID"].str.strip().eq("").any():
        raise ValueError("Missing ASV identifiers detected")

    if not counts["ASV_ID"].is_unique:
        raise ValueError("Duplicated ASV identifiers detected")

    # 2. Validate counts
    abundance = counts[sample_columns]

    if abundance.isna().any().any():
        raise ValueError("Missing values detected in ASV count table")

    numeric = abundance.apply(
        lambda column: column.map(parse_asv_count)
    )

    counts[sample_columns] = numeric.astype("int64")

    return counts


def find_zero_count_samples(counts: pd.DataFrame) -> list[str]:
    """Identify samples containing zero total ASV reads."""

    sample_columns = counts.columns.drop("ASV_ID")

    sample_totals = counts[sample_columns].sum(axis=0)

    return sample_totals[sample_totals == 0].index.tolist()

def load_dada2_taxonomy(path: str | Path) -> pd.DataFrame:
    """Load a DADA2 taxonomy table and validate ASV identifiers."""

    path = Path(path)

    validate_tsv_headers(path)

    taxonomy = pd.read_csv(path, sep="\t", dtype=str)

    if taxonomy.empty:
        raise ValueError("Taxonomy table is empty")

    if "ASV_ID" not in taxonomy.columns:
        raise ValueError("Missing required column: ASV_ID")

    if taxonomy.columns.duplicated().any():
        raise ValueError("Duplicated taxonomy columns detected")

    if taxonomy["ASV_ID"].isna().any():
        raise ValueError("Missing taxonomy ASV identifiers detected")

    if taxonomy["ASV_ID"].str.strip().eq("").any():
        raise ValueError("Empty taxonomy ASV identifiers detected")

    if not taxonomy["ASV_ID"].is_unique:
        raise ValueError("Duplicated taxonomy ASV identifiers detected")

    return taxonomy

def validate_asv_alignment(
    counts: pd.DataFrame,
    taxonomy: pd.DataFrame,
) -> None:
    """Verify that counts and taxonomy contain identical ASV identifiers."""

    count_ids = set(counts["ASV_ID"])
    taxonomy_ids = set(taxonomy["ASV_ID"])

    missing_taxonomy = count_ids - taxonomy_ids
    extra_taxonomy = taxonomy_ids - count_ids

    if missing_taxonomy:
        raise ValueError(
            f"Missing taxonomy for {len(missing_taxonomy)} ASVs"
        )

    if extra_taxonomy:
        raise ValueError(
            f"Taxonomy contains {len(extra_taxonomy)} unmatched ASVs"
        )