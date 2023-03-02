#' Peripheral Blood Mononuclear Cells (PBMCs)
#'
#' UMI counts for a subset of cells freely available from 10X Genomics
#'
#' @format A list with element 'counts' as sparse matrix (dgCMatrix, see Matrix package) of molecule counts,
#' and a meta_data data.frame with 'celltype' factor.
#' There are 600 rows (genes) and 300 columns (cells). This is a downsampled
#' version of a 10k healthy human PBMC dataset available from 10x Genomics.
#'
#' @source \url{https://www.10xgenomics.com/resources/datasets/10k-human-pbmcs-3-v3-1-chromium-controller-3-1-high}
"pbmc"
