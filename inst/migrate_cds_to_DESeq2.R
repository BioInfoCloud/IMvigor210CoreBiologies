#!/usr/bin/env Rscript
# =====================================================================
# migrate_cds_to_DESeq2.R  (DESeq-INDEPENDENT version)
# ---------------------------------------------------------------------
# Converts the legacy CountDataSet `cds` (old Bioconductor *DESeq* package,
# removed from Bioconductor >= 3.13 and no longer installable) into a
# modern DESeqDataSet, so that IMvigor210CoreBiologies no longer requires
# the uninstallable DESeq package.
#
# WHY THIS WORKS WITHOUT DESeq:
#   The package's own R/ functions only ever call DESeq2 / edgeR / limma.
#   The ONLY thing that depended on the old DESeq package was the shipped
#   data object `cds`, serialized as a CountDataSet. Even after removing
#   DESeq from DESCRIPTION, `data(cds)` would load an object whose class
#   is undefined. However, the CountDataSet slots (assayData, phenoData,
#   featureData) are reachable through S4 slot access without the DESeq
#   class definition being loaded, so we can extract counts + annotations
#   and rebuild a modern DESeqDataSet. No dead package required.
#
# PREREQUISITES (all modern, installable packages):
#   BiocManager::install(c("DESeq2", "SummarizedExperiment",
#                          "S4Vectors", "Biobase"))
#
# RUN FROM THE PACKAGE ROOT DIRECTORY:
#   Rscript inst/migrate_cds_to_DESeq2.R
#
# NOTE: the original data/cds.RData is backed up to data/cds.RData.orig
#       before it is overwritten.
# =====================================================================

pkg_root    <- if (requireNamespace("here", quietly = TRUE)) here::here() else getwd()
data_dir    <- file.path(pkg_root, "data")
orig_file   <- file.path(data_dir, "cds.RData")
backup_file <- file.path(data_dir, "cds.RData.orig")

## 0. sanity checks ----------------------------------------------------
stopifnot(file.exists(orig_file))

if (!requireNamespace("DESeq2", quietly = TRUE))
  stop("Please install DESeq2 first: BiocManager::install('DESeq2')")
if (!requireNamespace("SummarizedExperiment", quietly = TRUE))
  stop("Please install SummarizedExperiment: BiocManager::install('SummarizedExperiment')")
if (!requireNamespace("Biobase", quietly = TRUE))
  stop("Please install Biobase: BiocManager::install('Biobase')")

suppressMessages(library(DESeq2))
suppressMessages(library(Biobase))

## 1. backup -----------------------------------------------------------
if (!file.exists(backup_file)) {
  file.copy(orig_file, backup_file)
  message("Backed up original object to: ", backup_file)
}

## 2. load original CountDataSet (DESeq class is undefined -> benign warning) -
##    The wrapper form (load inside a function) would bind `cds` to the
##    wrong environment, so we load at top level and muffle the warning.
suppressWarnings(load(orig_file))
stopifnot(exists("cds"))

## 3. extract the three components via S4 slots (NO DESeq needed) -------
ad  <- cds@assayData
counts_mat <- if (is.environment(ad)) get("counts", envir = ad) else ad[["counts"]]
pd  <- as.data.frame(cds@phenoData@data)
fd  <- as.data.frame(cds@featureData@data)

stopifnot(is.matrix(counts_mat), nrow(counts_mat) > 0, ncol(counts_mat) > 0)

## 4. name-consistency sanity ------------------------------------------
if (!identical(rownames(counts_mat), rownames(fd)))
  stop("Feature (gene) names in counts and featureData do not align.")
if (!identical(colnames(counts_mat), rownames(pd)))
  stop("Sample names in counts and phenoData do not align.")

## 5. build a modern DESeqDataSet --------------------------------------
##    colData carries the former pData; rowData/mcols carry the former fData.
dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData   = S4Vectors::DataFrame(pd),
  design    = ~ 1
)
mcols(dds) <- S4Vectors::DataFrame(fd)

## best-effort carryover of legacy dispersion table (informational only)
tryCatch({
  dt <- cds@dispTable
  if (!is.null(dt) && length(dt) > 0) metadata(dds)$dispTable <- dt
}, warning = function(w) NULL, error = function(e) NULL)

## 6. overwrite the lazy-loaded dataset (must keep the name `cds`) ------
##    NOTE: `save(cds = dds, file = ...)` stores the object under the name
##    `dds` (save() uses the *argument names* as object names). We must
##    rebind to `cds` and save that symbol so `data(cds)` finds it.
cds <- dds
save(cds, file = orig_file)
message("Migrated `cds` -> DESeqDataSet and wrote: ", orig_file)
message("Dimensions: ", nrow(dds), " features x ", ncol(dds), " samples")
message("Access with counts(cds); colData(cds) [was pData]; rowData(cds) [was fData].")
message("NOTE: legacy code calling pData(cds)/fData(cds) should source inst/compatibility_shim.R")
