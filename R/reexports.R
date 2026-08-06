#' Re-exported accessors for the modern \code{DESeqDataSet} \code{cds}
#'
#' The shipped dataset \code{cds} was migrated from the legacy \code{CountDataSet}
#' (old Bioconductor \pkg{DESeq}, removed in Bioconductor >= 3.13) to a modern
#' \code{DESeqDataSet}. To let users call \code{counts(cds)}, \code{colData(cds)},
#' etc. immediately after \code{library(IMvigor210CoreBiologies)} -- without also
#' having to attach \pkg{DESeq2}/\pkg{SummarizedExperiment} -- these accessors are
#' re-exported here.
#' @name reexports
#' @keywords internal
NULL

#' @importFrom DESeq2 counts
#' @export
counts <- DESeq2::counts

#' @importFrom SummarizedExperiment colData
#' @export
colData <- SummarizedExperiment::colData

#' @importFrom SummarizedExperiment rowData
#' @export
rowData <- SummarizedExperiment::rowData

#' @importFrom SummarizedExperiment assay
#' @export
assay <- SummarizedExperiment::assay

#' @importFrom SummarizedExperiment assays
#' @export
assays <- SummarizedExperiment::assays

#' @importFrom S4Vectors metadata
#' @export
metadata <- S4Vectors::metadata
