# ASV Explorer

ASV Explorer is a Python and Nextflow project for validating,
exploring, and visualizing taxonomic abundance data.

It currently supports two complementary workflows:

- Taxonomic abundance visualization using Python, including
  nested donut charts with collision-reduced labels.
- DADA2 ASV count validation using Nextflow DSL2 and Docker,
  including SILVA taxonomy alignment and quality-control reports.

The visualization example uses aggregated transcript abundance
(TPM) from the INTERES WP1 dataset. The repository does not
contain the original 77 MB processed dataset; it includes only
the aggregated values needed to reproduce the example figure.

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"

asv-explorer
pytest
```

Label smaller taxa in the outer ring by lowering the percentage threshold:

```bash
asv-explorer --min-label-percent 0.5
```

The command creates:

- `results/taxonomic_abundance_summary.csv`
- `results/taxonomic_abundance_donut.png`

## ASV validation with Nextflow

ASV Explorer v0.2.0 includes a Nextflow DSL2 workflow for validating
DADA2 ASV count tables and their corresponding taxonomy assignments.

### Requirements

- Nextflow (tested with 26.10.4)
- Java 21
- Docker

### Build the container

```bash
docker build -f docker/Dockerfile -t asv-explorer:0.2.0 .
```

### Run the workflow

```bash
nextflow run main.nf \
  --counts path/to/ASV_table.tsv \
  --taxonomy path/to/ASV_taxonomy.tsv \
  --outdir results/validation \
  -with-docker
```

The workflow generates:

- `validation.json`: dataset dimensions, read conservation,
  taxonomy alignment, and input provenance.
- `sample_summary.tsv`: total reads and observed ASVs per sample.
- `feature_summary.tsv`: total abundance and sample prevalence
  per ASV.

### Input requirements

The count table must contain an `ASV_ID` column followed by
sample columns containing non-negative integer read counts.

The taxonomy table must contain an `ASV_ID` column.

Both tables must contain matching, unique ASV identifiers.

Missing taxonomic classifications are permitted, but missing
ASV identifiers or unmatched ASVs are rejected.

Zero-count samples are reported without being removed.

The validation workflow does not perform taxonomic filtering,
normalization, or ecological diversity analysis.

### Testing

Python unit tests:

```bash
python -m pytest -q --ignore=tests/test_nextflow.py
```

Nextflow integration test (requires the Docker image):

```bash
python -m pytest tests/test_nextflow.py -v
```

## Input format

The default CSV schema is intentionally generic:

```csv
parent_taxon,taxon,abundance
Dinoflagellata,Gonyaulacales,1250.4
Ochrophyta,Pseudo-nitzschia,830.2
```

Repeated taxonomic groups are summed automatically. Abundance values must be numeric
and non-negative.

Column names can be adapted without modifying the code:

```bash
asv-explorer \
  --input path/to/table.csv \
  --parent-column Class \
  --taxon-column Final_Taxonomy \
  --abundance-column tpm \
  --output-dir results
```

## Roadmap: ASV analysis

Completed in v0.2.0:

- DADA2 ASV count table import and validation
- SILVA taxonomy table import and ASV alignment
- Detection of zero-count samples
- Read conservation checks
- Sample and feature summaries
- Nextflow DSL2 execution with Docker

Planned:

- Per-sample relative abundance
- Configurable ASV filtering with read accounting
- Taxonomic aggregation by rank
- Stacked taxonomic abundance plots
- Alpha and beta diversity analysis
- Support for additional input formats and taxonomic databases

## Data provenance

The example CSV is a taxonomic aggregation of these three source columns:
`Class`, `Final_Taxonomy`, and `tpm`. See `scripts/prepare_example_data.R` for the
reproducible preparation step. The full local source remains outside Git.

## License

MIT
