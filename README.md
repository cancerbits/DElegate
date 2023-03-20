# DElegate R package

## Wrapper and helper functions to use bulk RNA-seq differential expression methods with single-cell data

Developed by Christoph Hafemeister in the [Developmental Cancer Genomics](https://cancerbits.github.io) group at [St. Anna Children's Cancer Research Institute (CCRI)](https://ccri.at)

## Quick start

Install from GitHub
```r
remotes::install_github('cancerbits/DElegate')
```

Given a Seurat object `s`, run differential expression tests between each cluster and the rest of the cells.
```r
de_results <- DElegate::findDE(object = s)
```

Or find all cluster markers and show top 5 for each cluster
```r
marker_results <- DElegate::FindAllMarkers2(object = s)
dplyr::filter(marker_results, feature_rank < 6)
```

## Vignette

An overview of the functionality, including examples, can be found [here](https://cancerbits.github.io/projects/introduction_to_DElegate.html)

## Summary

`DElegate` is an R package that allows bulk RNA-seq differential expression methods to be used with single-cell data. It is a light wrapper around 
[`DESeq2`](https://doi.org/doi:10.18129/B9.bioc.DESeq2), 
[`edgeR`](https://doi.org/doi:10.18129/B9.bioc.edgeR), and 
[`limma`](https://doi.org/doi:10.18129/B9.bioc.limma), similar to the [`Libra`](https://github.com/neurorestore/Libra) package. In contrast to `Libra`, `DElegate` focuses on a few DE methods and will assign cells to pseudo-replicates if no true replicates are available.


All `DElegate` functionality is contained in one function - `findDE()`. It has one mandatory input argument: `object`, which may be of class

* `Seurat` - the count matrix will be extracted from the `'RNA'` assay
* `SingleCellExperiment` - the count matrix will be extracted via `counts()`
* `dgCMatrix` - sparse matrix of the `Matrix` package
* `matrix`

To indicate the cell group memberships, you have several options, depending on input type:

* `Seurat` - in the object via `Idents(object)`, or use the `group_column` argument
* `SingleCellExperiment` - in the object via `colLables(object)`, or use the `group_column` argument
* `dgCMatrix`, or `matrix` - use the `meta_data` and `group_column` arguments

`DElegate` uses bulk RNA-seq DE methods and relies on replicates. If no true replicates are available, it assigns cells to pseudo-replicates. However, if replicates are available in the input, the `replicate_column` argument can be used to indicate where to find the labels.

To tell `findDE()` which cell groups to compare, use the `compare` argument. We provide several ways to set up the comparisons that will be tested:

* `'each_vs_rest'`, the default, does multiple comparisons, one per group vs all remaining cells
* `'all_vs_all'`, also does multiple comparisons, covering all group pairs
* a length one character vector, e.g. `'MONOCYTES'`, does one comparison between that group and the remaining cells
* a length two character vector, e.g. `c('T CELLS', 'B CELLS')`, does one comparison between those two groups
* a list of length two, e.g. `list(c('T CELLS', 'B CELLS'), c('MONOCYTES'))`, does one comparison after combining groups

Finally, there are currently three DE methods supported

* `'edger'` uses `edgeR::glmQLFit`
* `'deseq'` uses `DESeq2::DESeq(test = 'Wald')`
* `'limma'` uses `limma::eBayes(trend = TRUE, robust = TRUE)`

For complete details, consult the package documentation: `?DElegate::findDE`. 

## Parallelization

`DElegate` supports parallelization via the `future` package. 
For example, to use the multicore strategy with 12 workers you may set
`future::plan(strategy = 'future::multicore', workers = 12)`.
See more details at the [`future` website](https://future.futureverse.org)

Note that every comparison is run single-threaded, but multiple comparisons will be done in parallel.

## Progress updates

For reporting progress updates, `DElegate` relies on the `progressr` package. By default no progress updates are rendered, but may be turned on in an R session: `progressr::handlers(global = TRUE)` and its default presentation modified (e.g. `progressr::handlers(progressr::handler_progress)`). 
See more details at the [`progressr` website](https://progressr.futureverse.org)

