#!/usr/bin/env Rscript

# Create the small, publishable example used by the Python application.
# This script is only needed to regenerate the example from the local source.

args <- commandArgs(trailingOnly = TRUE)
source_path <- if (length(args) >= 1) args[[1]] else {
  "/mnt/data/ProjectsData/bacterial-suppression/data/processed/metatranscriptomics/INTERES_Euk_WP1_TPMs_Annotated_1KEGG_ko.parquet"
}
output_path <- if (length(args) >= 2) args[[2]] else {
  "data/example_taxonomic_abundance.csv"
}

required <- c("Class", "Final_Taxonomy", "tpm")
table <- arrow::read_parquet(source_path, col_select = required, as_data_frame = TRUE)

table$Class[is.na(table$Class) | table$Class == ""] <- "Unclassified"
table$Final_Taxonomy[is.na(table$Final_Taxonomy) | table$Final_Taxonomy == ""] <- "Unclassified"

example <- aggregate(
  table$tpm,
  by = list(parent_taxon = table$Class, taxon = table$Final_Taxonomy),
  FUN = sum,
  na.rm = TRUE
)
names(example)[3] <- "abundance"
example <- example[example$abundance > 0, ]
example <- example[order(example$abundance, decreasing = TRUE), ]

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write.csv(example, output_path, row.names = FALSE, quote = TRUE)
cat("Created", output_path, "with", nrow(example), "taxonomic groups\n")
