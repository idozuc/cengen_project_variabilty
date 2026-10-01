#!/usr/bin/env Rscript

# Create gene and neuron classification tables for enrichment analysis.
script_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_file[[1L]]))
project_root <- dirname(dirname(dirname(script_path)))
source(file.path(project_root, "src", "R", "project_io.R"))
options <- parse_named_arguments(
  commandArgs(trailingOnly = TRUE),
  list(master = NA_character_, wormcat = NA_character_, cengen_supplement = NA_character_, output_dir = NA_character_),
  paste(
    "Usage: Rscript scripts/analysis/01_prepare_mappings.R --master FILE.csv.gz",
    "--wormcat FILE.csv --cengen-supplement FILE.xlsx --output-dir PATH [--config PATH]"
  )
)
master_file <- normalizePath(require_path_option(options$master, "--master"), mustWork = TRUE)
wormcat_file <- normalizePath(require_path_option(options$wormcat, "--wormcat"), mustWork = TRUE)
cengen_file <- normalizePath(require_path_option(options$cengen_supplement, "--cengen-supplement"), mustWork = TRUE)
output_dir <- require_path_option(options$output_dir, "--output-dir", must_exist = FALSE)

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop('Install readxl with install.packages("readxl")')
}

# Stop when a source table lacks columns required below.
check_columns <- function(x, required, path) {
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop("Missing columns in ", path, ": ", paste(missing, collapse = ", "))
  }
}

master <- read.csv(gzfile(master_file), stringsAsFactors = FALSE)
check_columns(master, c("stable_id", "gene_name", "cell_type"), master_file)

# Match each BASiCS gene to its WormCat 2.0 hierarchy.
wormcat <- read.csv(wormcat_file, stringsAsFactors = FALSE, check.names = FALSE)
wormcat_columns <- c(
  "Wormbase ID", "Category 1", "Category 2", "Category 3",
  "Automated Description"
)
check_columns(wormcat, wormcat_columns, wormcat_file)
if (anyDuplicated(wormcat[["Wormbase ID"]])) {
  stop("WormCat contains duplicate WormBase IDs")
}

genes <- unique(master[c("stable_id", "gene_name")])
index <- match(genes$stable_id, wormcat[["Wormbase ID"]])
gene_classes <- data.frame(
  stable_id = genes$stable_id,
  gene_name = genes$gene_name,
  category_1 = wormcat[["Category 1"]][index],
  category_2 = wormcat[["Category 2"]][index],
  category_3 = wormcat[["Category 3"]][index],
  description = wormcat[["Automated Description"]][index],
  stringsAsFactors = FALSE
)
gene_classes <- gene_classes[!is.na(index), ]
gene_classes <- gene_classes[order(gene_classes$stable_id), ]
rownames(gene_classes) <- NULL

# Extract direct neuron-family mappings from CeNGen Supplement 13B-E.
neuron_sheets <- readxl::excel_sheets(cengen_file)[2:5]
read_neurons <- function(sheet) {
  x <- readxl::read_excel(
    cengen_file,
    sheet = sheet,
    skip = 2,
    col_names = FALSE,
    col_types = "text",
    .name_repair = "minimal"
  )
  x <- as.data.frame(x[, 1:5])
  names(x) <- c("gene", "cell_type", "tpm", "proportion", "neuron_family")
  x[c("cell_type", "neuron_family")]
}
neuron_classes <- unique(do.call(rbind, lapply(neuron_sheets, read_neurons)))
neuron_classes <- neuron_classes[
  !is.na(neuron_classes$cell_type) & !is.na(neuron_classes$neuron_family),
]

counts <- table(neuron_classes$cell_type)
if (any(counts > 1L)) stop("A neuron has more than one published family")

# Keep the complete published neuron-to-functional-group lookup.
all_neuron_classes <- neuron_classes[order(neuron_classes$cell_type), ]
names(all_neuron_classes) <- c("neuron", "neuron_functional_group")
rownames(all_neuron_classes) <- NULL

cell_types <- unique(master$cell_type)
neuron_classes <- neuron_classes[neuron_classes$cell_type %in% cell_types, ]
neuron_classes <- neuron_classes[order(neuron_classes$cell_type), ]
rownames(neuron_classes) <- NULL

# Save only mappings supported by the selected sources.
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(
  gene_classes,
  file.path(output_dir, "gene_classes.csv"),
  row.names = FALSE,
  na = "NA"
)
write.csv(
  neuron_classes,
  file.path(output_dir, "neuron_classes.csv"),
  row.names = FALSE,
  na = "NA"
)
write.csv(
  all_neuron_classes,
  file.path(output_dir, "neuron_functional_groups.csv"),
  row.names = FALSE,
  na = "NA"
)

cat("Mapped", nrow(gene_classes), "of", nrow(genes), "genes with WormCat 2.0\n")
cat(
  "Mapped", nrow(neuron_classes), "of", length(cell_types),
  "analyzed cell types with CeNGen\n"
)
cat("Saved", nrow(all_neuron_classes), "published neuron mappings\n")
