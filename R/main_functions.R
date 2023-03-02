
#' Differential expression tests for single-cell RNA-seq data
#'
#' @param object Input object (of class Seurat, SingleCellExperiment, dgCMatrix, or matrix)
#' @param meta_data Data frame of meta data; can be NULL if object is Seurat or
#' SingleCellExperiment and cell identities are set; default is NULL
#' @param group_column String or integer indicating what column of meta_data to use as group indicator;
#' if NULL, the object cell identities are used; default is NULL
#' @param replicate_column String or integer indicating what column of meta_data
#' to use as replicate indicator; default is NULL
#' @param compare Specifies which groups to compare - see details; default is 'each_vs_rest'
#' @param method DE method to use; default is edger
#' @param verbosity Integer controlling how many messages the function prints;
#' 0 is silent, 1 prints some messages, 2 prints more messages
#' @returns A data frame of results
#'
#' @section Details:
#' Compare groups of cells using DESeq2, edgeR, or limma-trend.
#'
#' There are multiple ways the group comparisons can be specified based on the compare
#' parameter. The default, \code{'each_vs_rest'}, does multiple comparisons, one per
#' group vs all remaining cells. \code{'all_vs_all'}, also does multiple comparisons,
#' covering all group pairs. If compare is set to a length two character vector, e.g.
#' \code{c('T-cells', 'B-cells')}, one comparison between those two groups is done.
#' To put multiple groups on either side of a single comparison, use a list of length two.
#' E.g. \code{compare = list(c('cluster1', 'cluster5'), c('cluster3'))}.
#'
#' @export
#'
findDE <- function(object,
                   meta_data = NULL,
                   group_column = NULL,
                   replicate_column = NULL,
                   compare = 'each_vs_rest',
                   method = 'edger',
                   verbosity = 1) {
  # extract the data from the input object
  de_data <- get_data(object, meta_data, group_column, replicate_column, verbosity)
  # set up the comparisons
  group_levels = levels(x = de_data$grouping)
  comparisons <- set_up_comparisons(group_levels = group_levels, compare = compare, verbosity = verbosity)
  print_comparisons(comparisons, verbosity)
  # run DE
  run_de_comparisons(counts = de_data$counts,
                     grouping = de_data$grouping,
                     replicate_label = de_data$replicate_label,
                     comparisons = comparisons,
                     method = method,
                     verbosity = verbosity)
}

#' Find markers for all groups
#'
#' This calls [findDE()] with compare = 'each_vs_rest' and filters and orders the results.
#'
#' @inheritParams findDE
#' @returns A data frame of results
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @export
#'
FindAllMarkers2 <- function(object,
                            meta_data = NULL,
                            group_column = NULL,
                            replicate_column = NULL,
                            method = 'edger',
                            verbosity = 1) {
  res <- findDE(object = object,
                meta_data = meta_data,
                group_column = group_column,
                replicate_column = replicate_column,
                compare = 'each_vs_rest',
                method = method,
                verbosity = verbosity)
  if (method == 'edger') {
    res <- dplyr::filter(res, .data$logFC < 0) %>%
      dplyr::group_by(.data$group1) %>%
      dplyr::arrange(.data$group1, .data$PValue, -.data$F, .by_group = TRUE) %>%
      dplyr::mutate(feature_rank = 1:dplyr::n()) %>%
      dplyr::select(-.data$group2)
  }
  if (method == 'deseq') {
    res <- dplyr::filter(res, !is.na(.data$log2FoldChange), .data$log2FoldChange < 0) %>%
      dplyr::group_by(.data$group1) %>%
      dplyr::arrange(.data$group1, .data$pvalue, .data$stat, .by_group = TRUE) %>%
      dplyr::mutate(feature_rank = 1:dplyr::n()) %>%
      dplyr::select(-.data$group2)
  }
  if (method == 'limma') {
    res <- dplyr::filter(res, .data$logFC < 0) %>%
      dplyr::group_by(.data$group1) %>%
      dplyr::arrange(.data$group1, .data$P.Value, -.data$B, .by_group = TRUE) %>%
      dplyr::mutate(feature_rank = 1:dplyr::n()) %>%
      dplyr::select(-.data$group2)
  }

  return(as.data.frame(res))
}

