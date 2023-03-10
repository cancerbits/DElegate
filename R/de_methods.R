# DE functions are defined here
# each one expects a count matrix (rows are genes)
# and a grouping factor (two levels assumed)

#' @importFrom magrittr %>%
run_de_simple <- function(counts, grouping, replicate_label, method, order_results, verbosity) {
  if (!any(method %in% c('deseq', 'edger', 'limma')) | length(method) != 1) {
    stop("method needs to be of length one and in c('deseq', 'edger', 'limma')")
  }

  # we are only dealing with bulk methods, so create pseudobulk data
  pb <- make_pb(counts, grouping, replicate_label)
  if (verbosity > 1) {
    message(paste0(c('group vs replicate table', utils::capture.output(pb$md)), collapse = "\n"))
  }

  # the DE methods we use, treat the first group level as reference
  # we want this to be group2, so relevel here
  pb$md$grouping <- stats::relevel(pb$md$grouping, ref = levels(pb$md$grouping)[2])

  if (method == 'deseq') {
    if (!requireNamespace("DESeq2", quietly = TRUE)) {
      stop('DESeq2 package not found - please install it to use "method = deseq"')
    }
    res <- run_deseq_simple(pb$counts, pb$md$grouping, order_results)
  }
  if (method == 'edger') {
    if (!requireNamespace("edgeR", quietly = TRUE)) {
      stop('edgeR package not found - please install it to use "method = edger"')
    }
    res <- run_edger_simple(pb$counts, pb$md$grouping, order_results)
  }
  if (method == 'limma') {
    if (!requireNamespace("edgeR", quietly = TRUE) | !requireNamespace("limma", quietly = TRUE)) {
      stop('edgeR or limma package not found - please install both to use "method = limma"')
    }
    res <- run_limma_simple(pb$counts, pb$md$grouping, order_results)
  }
  return(res)
}


run_deseq_simple <- function(counts, grouping, order_results) {
  dds <- DESeq2::DESeqDataSetFromMatrix(counts, data.frame(grouping = grouping), ~ grouping)
  dds <- DESeq2::DESeq(dds, test = 'Wald', quiet = TRUE)
  res <- DESeq2::results(dds, tidy = TRUE) %>%
    dplyr::rename('feature' = 'row')
  if (order_results) {
    res <- dplyr::arrange(res, .data$pvalue, -abs(.data$stat))
  }
  return(res)
}

#' @importFrom magrittr %>%
#' @importFrom rlang .data
run_edger_simple <- function(counts, grouping, order_results) {
  design <- stats::model.matrix(~ grouping)

  y <- edgeR::DGEList(counts = counts, group = grouping)
  y <- edgeR::calcNormFactors(y)
  y <- edgeR::estimateDisp(y, design)
  fit <- edgeR::glmQLFit(y, design)
  test_res <- edgeR::glmQLFTest(fit)

  res <- test_res$table %>%
    tibble::rownames_to_column(var = 'feature') %>%
    dplyr::mutate(AdjPValue = stats::p.adjust(.data$PValue, method = 'fdr'))
  if (order_results) {
    res <- dplyr::arrange(res, .data$PValue, -.data$F)
  }
  return(res)
}

run_limma_simple <- function(counts, grouping, order_results) {
  design <- stats::model.matrix(~ grouping)
  n <- ncol(design)
  allow_trend <- TRUE
  y <- edgeR::DGEList(counts = counts, group = grouping)
  y <- edgeR::calcNormFactors(y)
  y <- edgeR::cpm(y, log = TRUE, prior.count = 3)
  fit <- limma::lmFit(y, design)
  fit <- limma::eBayes(fit, trend = allow_trend, robust = allow_trend)
  res <- data.frame(feature = rownames(fit$p.value),
                    logFC = fit$coefficients[, n],
                    AveExpr = fit$Amean,
                    t = fit$t[, n],
                    P.Value = fit$p.value[, n],
                    adj.P.Val = stats::p.adjust(fit$p.value[, n], method = 'fdr'),
                    B = fit$lods[, n],
                    row.names = NULL)
  if (order_results) {
    res <- dplyr::arrange(res, .data$P.Value, -.data$B)
  }
  return(res)
}
