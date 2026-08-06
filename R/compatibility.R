#' Register legacy \code{pData}/\code{fData} methods for \code{DESeqDataSet}
#'
#' Historically \code{cds} was a \code{CountDataSet} (Bioconductor \pkg{DESeq})
#' and analysis code called \code{pData(cds)} / \code{fData(cds)}. After the
#' migration to \code{DESeqDataSet} those accessors are \code{colData} /
#' \code{rowData}. We register drop-in \code{pData}/\code{fData} methods on
#' package load so the legacy API keeps working without any extra steps.
#' (A standalone copy of this logic is also available as
#' \code{inst/compatibility_shim.R}.)
.onLoad <- function(libname, pkgname) {
  if (requireNamespace("Biobase", quietly = TRUE) &&
      requireNamespace("DESeq2", quietly = TRUE)) {
    methods::setMethod(
      "pData", "DESeqDataSet",
      function(object) as.data.frame(colData(object))
    )
    methods::setMethod(
      "fData", "DESeqDataSet",
      function(object) as.data.frame(rowData(object))
    )
  }
}
