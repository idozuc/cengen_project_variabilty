#!/usr/bin/env Rscript

# Fraction of cells with >= 1 UMI for selected genes in each cluster of the
# accepted HSN, SMD and NSM partitions (Figure 12).
#
# Writes OUTPUT_DIR/data/subgroup_markers.csv

script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_file[[1L]]))
project_root <- dirname(dirname(dirname(script_path)))
source(file.path(project_root, "src", "R", "project_io.R"))
source(file.path(project_root, "src", "R", "figures.R"))
source(file.path(project_root, "src", "R", "clustering.R"))

options <- figure_options(
  "Usage: Rscript scripts/figures/03_export_subgroup_markers.R --config config/figures.yaml"
)

genes <- list(
  HSN = c("unc-86", "ast-1", "hlh-3", "egl-46", "tph-1", "cat-1", "bas-1", "flp-19", "nlp-3", "rab-3"),
  SMD = c("flp-5", "tbx-2", "pdf-1", "nlp-21", "flp-7", "egl-21", "sbt-1", "rab-3"),
  NSM = c("del-3", "del-7", "tts-2", "pgp-14", "ceh-2", "tph-1", "cat-1", "mod-5", "rab-3")
)

sce <- read_sce(options)
summary_table <- read.csv(input_file(options$clustering_dir, "run_summary.csv"))

rows <- list()
for (cell_type in names(genes)) {
  # Same cells as the clustering run: experiments with >= 5 cells.
  selected <- filter_rare_experiments(select_cell_type(sce, cell_type), min_cells = 5L)$sce
  g <- summary_table$primary_g[summary_table$cell_type == cell_type]
  assignments <- read.csv(input_file(options$clustering_dir, cell_type, "assignments_by_resolution.csv"))
  cluster <- setNames(
    assignments[[paste0("cluster_g", g)]],
    assignments$observation
  )[colnames(selected)]
  symbol <- as.character(rowData(selected)$gene_short_name)

  for (gene in genes[[cell_type]]) {
    detected <- counts(selected)[which(symbol == gene)[1L], ] > 0
    for (k in sort(unique(cluster))) {
      rows[[length(rows) + 1L]] <- data.frame(
        cell_type = cell_type,
        gene = gene,
        cluster = k,
        cells = sum(cluster == k),
        detection_fraction = mean(detected[cluster == k])
      )
    }
  }
}

out <- do.call(rbind, rows)
write.csv(out, file.path(options$data_dir, "subgroup_markers.csv"), row.names = FALSE)
print(out)
