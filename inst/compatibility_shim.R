# =====================================================================
# compatibility_shim.R
# ---------------------------------------------------------------------
# After migrating `cds` to a DESeqDataSet, legacy analysis code that calls
# pData(cds) / fData(cds) would break, because those Biobase accessors are
# defined for eSet objects, not DESeqDataSet. This shim adds drop-in
# methods so the original API keeps working:
#
#   counts(cds)  -> unchanged (DESeqDataSet already supports it)
#   pData(cds)   -> colData(cds)
#   fData(cds)   -> rowData(cds)
#
# Source it once in your session after `library(IMvigor210CoreBiologies)`:
#   source(system.file("compatibility_shim.R", package = "IMvigor210CoreBiologies"))
# =====================================================================

suppressMessages(library(DESeq2))
suppressMessages(library(Biobase))

## Biobase defines the `pData`/`fData` generics with signature (object),
## so the methods MUST match that signature exactly (no `...`).
setMethod("pData", "DESeqDataSet",
          function(object) as.data.frame(colData(object)))

setMethod("fData", "DESeqDataSet",
          function(object) as.data.frame(rowData(object)))
