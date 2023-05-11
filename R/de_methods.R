# DE functions are defined here
# each one expects a count matrix (rows are genes)
# and a grouping factor (two levels assumed)

#' @importFrom magrittr %>%
run_de_simple <- function(counts, grouping, replicate_label, method, order_results, lfc_shrinkage, verbosity) {
  if (!any(method %in% c('deseq', 'edger', 'limma')) | length(method) != 1) {
    stop("method needs to be of length one and in c('deseq', 'edger', 'limma')")
  }

  # we are only dealing with bulk methods, so create pseudobulk data
  pb <- make_pb(counts, grouping, replicate_label)
  if (verbosity > 1) {
    message(paste0(c('group vs replicate table', utils::capture.output(pb$md)), collapse = "\n"))
  }

  if (verbosity > 0 && !is.null(lfc_shrinkage) && method == 'deseq') {
    message(sprintf('calculate shrunken log2 fold changes of type %s', lfc_shrinkage))
  }

  # the DE methods we use, treat the first group level as reference
  # we want this to be group2, so relevel here
  pb$md$grouping <- stats::relevel(pb$md$grouping, ref = levels(pb$md$grouping)[2])

  if (method == 'deseq') {
    if (!requireNamespace("DESeq2", quietly = TRUE)) {
      stop('DESeq2 package not found - please install it to use "method = deseq"')
    }
    res <- run_deseq_simple(pb$counts, pb$md$grouping, order_results, lfc_shrinkage)
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

  # add some feature-level stats
  det_rate <- detection_rate(counts, grouping)[res$feature, ]
  colnames(det_rate) <- c('rate1', 'rate2')
  rownames(det_rate) <- NULL
  res <- cbind(res, det_rate)

  return(res)
}


run_deseq_simple <- function(counts, grouping, order_results, lfc_shrinkage) {
  dds <- DESeq2::DESeqDataSetFromMatrix(counts, data.frame(grouping = grouping), ~ grouping)
  dds <- DESeq2::DESeq(dds, test = 'Wald', quiet = TRUE)
  if (is.null(lfc_shrinkage)) {
    res <- DESeq2::results(dds)
  }
  else {
    if (lfc_shrinkage %in% c("apeglm", "ashr", "normal")) {
      coef_name <- DESeq2::resultsNames(dds)[2]
      res <- DESeq2::lfcShrink(dds, coef = coef_name, type = lfc_shrinkage)
    }
    else {
      stop('lfc_shrinkage should be set to NULL or one of "apeglm", "ashr", "normal"')
    }
  }
  res <- as.data.frame(res) %>%
    tibble::rownames_to_column(var = 'feature') %>%
    dplyr::rename('ave_expr' = 'baseMean', 'log_fc' = 'log2FoldChange')
  # some shrinkage methods drop the stat column - add NA
  if (!'stats' %in% colnames(res)) {
    res <- tibble::add_column(res, stat = NA_real_, .before = 'pvalue')
  }
  # select the columns we want to return
  res <- dplyr::select(res, .data$feature, .data$ave_expr, .data$log_fc, .data$stat, .data$pvalue, .data$padj)
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
    dplyr::rename('ave_expr' = 'logCPM', 'log_fc' = 'logFC', 'stat' = 'F', 'pvalue' = 'PValue') %>%
    dplyr::mutate(padj = stats::p.adjust(.data$pvalue, method = 'fdr')) %>%
    dplyr::select(.data$feature, .data$ave_expr, .data$log_fc, .data$stat, .data$pvalue, .data$padj)
  if (order_results) {
    res <- dplyr::arrange(res, .data$pvalue, -.data$stat)
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
                    ave_expr = fit$Amean,
                    log_fc = fit$coefficients[, n],
                    stat = fit$lods[, n],
                    pvalue = fit$p.value[, n],
                    padj = stats::p.adjust(fit$p.value[, n], method = 'fdr'),
                    row.names = NULL)
  res <- tibble::as_tibble(res)
  if (order_results) {
    res <- dplyr::arrange(res, .data$pvalue, -.data$stat)
  }
  return(res)
}
