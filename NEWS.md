# DElegate 1.2.0

## New features

* `findDE()` gains a `compare_is_ref` argument that indicates that the `compare` parameter defines the reference group levels. If set, all non-reference groups will be compated to the reference group. Default is `FALSE` for same behaviour as in previous version.


# DElegate 1.1.0

## New features

* `findDE()` gains an `lfc_shrinkage` argument that controls which type of logFC-shrinkage is used with the *DESeq2* method. Default is `NULL` for no shrinkage as in previous version.


# DElegate 1.0.0

* First feature-complete release.
* Accompanies the pre-print [*Single-cell RNA-seq differential expression tests within a sample should use pseudo-bulk data of pseudo-replicates*](https://doi.org/10.1101/2023.03.28.534443).
