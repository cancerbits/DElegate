test_that("findDE works", {
  skip_on_cran()
  set.seed(42)

  # default with sparse matrix input
  counts <- DElegate::pbmc$counts
  meta_data <- DElegate::pbmc$meta_data
  res <- findDE(object = counts, meta_data = meta_data, group_column = 'celltype')
})
