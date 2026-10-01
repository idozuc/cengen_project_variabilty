source(file.path("src", "R", "project_io.R"))
source(file.path("src", "R", "stress_gene_sets.R"))

test_that("stress sets are stable and nested", {
  expect_length(stress_genes_short_14, 14L)
  expect_length(stress_genes_full_199, 199L)
  expect_true(all(stress_genes_short_14 %in% stress_genes_full_199))
})

test_that("BASiCS result directories use the canonical layout", {
  expect_identical(
    basics_result_dir("run", 2L, "SMD"),
    file.path("run", "cell_types", "002_SMD")
  )
})
