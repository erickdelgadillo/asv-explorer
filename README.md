# ASV Explorer

A small Python project for exploring and visualizing taxonomic relative abundance.
It turns a simple abundance table into a validated summary and a nested-donut chart:
the inner ring represents a higher taxonomic level and the outer ring represents the
more specific taxa.

The current example uses aggregated transcript abundance (TPM) from the INTERES WP1
dataset. The repository does **not** contain the original 77 MB processed dataset;
it includes only the aggregated values needed to reproduce the example figure.

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"

asv-explorer
pytest
```

The command creates:

- `results/taxonomic_abundance_summary.csv`
- `results/taxonomic_abundance_donut.png`

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

The next iteration will add an ASV-oriented loader supporting:

- ASV count/abundance matrices plus a taxonomy table
- per-sample relative abundance
- filtering of rare ASVs
- aggregation by phylum, class, order, family, or genus
- stacked bar charts alongside the current donut chart

## Data provenance

The example CSV is a taxonomic aggregation of these three source columns:
`Class`, `Final_Taxonomy`, and `tpm`. See `scripts/prepare_example_data.R` for the
reproducible preparation step. The full local source remains outside Git.

## License

MIT
