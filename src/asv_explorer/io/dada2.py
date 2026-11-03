"""Input and validation functions for DADA2 ASV count tables."""

from pathlib import Path

import numpy as np
import pandas as pd


def load_dada2_counts(path: str | Path) -> pd.DataFrame:
    """Load and validate a DADA2 ASV count table."""

    path = Path(path)

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

    try:
        numeric = abundance.apply(pd.to_numeric, errors="raise")
    except (ValueError, TypeError) as exc:
        raise ValueError("Non-numeric ASV counts detected") from exc

    values = numeric.to_numpy(dtype=float)

    if not np.isfinite(values).all():
        raise ValueError("Non-finite ASV counts detected")

    if (values < 0).any():
        raise ValueError("Negative ASV counts detected")

    if not np.equal(values, np.floor(values)).all():
        raise ValueError("Fractional ASV counts detected")

    # 3. Return validated integer counts
    try:
        numeric = numeric.astype("int64")
    except (ValueError, OverflowError) as exc:
        raise ValueError("ASV counts exceed supported integer range") from exc

    counts[sample_columns] = numeric

    return counts


def find_zero_count_samples(counts: pd.DataFrame) -> list[str]:
    """Identify samples containing zero total ASV reads."""

    sample_columns = counts.columns.drop("ASV_ID")

    sample_totals = counts[sample_columns].sum(axis=0)

    return sample_totals[sample_totals == 0].index.tolist()

def load_dada2_taxonomy(path: str | Path) -> pd.DataFrame:
    """Load a DADA2 taxonomy table and validate ASV identifiers."""

    path = Path(path)

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