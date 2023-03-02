test_that("findDE works", {
  skip_on_cran()
  set.seed(42)

  # default with sparse matrix input
  counts <- scDE::pbmc$counts
  meta_data <- scDE::pbmc$meta_data
  res <- findDE(object = counts, meta_data = meta_data, group_column = 'celltype')
})
