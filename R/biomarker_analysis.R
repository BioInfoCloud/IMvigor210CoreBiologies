#' Immunotherapy-biomarker helper functions for the IMvigor210 \code{cds}
#'
#' A small, self-contained toolkit that wraps the most common biomarker
#' analyses performed on the bundled \code{cds} (IMvigor210 / atezolizumab
#' metastatic urothelial cancer cohort). Every function takes the \code{cds}
#' object (and a gene symbol) and returns a \pkg{ggplot2} / \pkg{survminer}
#' plot object that can be printed, saved, or further customized. Expression is
#' computed as \code{log2(TPM + 1)} (consistent with the reference
#' \code{RNH1分析.R}) with an edgeR low-expression filter.
#'
#' @section Cohorts:
#' Three cohorts are used internally, mirroring the original analysis:
#' \describe{
#'   \item{ALL}{all 348 samples (incl. 50 response-unknown / NE); for analyses
#'     that do not need a response label (correlation, all-sample heatmap,
#'     expression-level KM).}
#'   \item{EVAL}{298 response-known samples (68 Responder / 230 Non-responder);
#'     for every response-dependent analysis.}
#'   \item{IMM}{244 samples with complete immune annotation (61 Responder /
#'     183 Non-responder); for immune-microenvironment context.}
#' }
#'
#' @importFrom pheatmap pheatmap
#' @importFrom survminer ggsurvplot surv_pvalue
#' @importFrom pROC roc auc ci.auc
#' @importFrom ggsignif geom_signif
#' @importFrom ggrepel geom_text_repel
#' @importFrom patchwork plot_layout
#' @importFrom edgeR DGEList calcNormFactors cpm
#' @noRd
NULL

##' Default immune-checkpoint / cytotoxic gene panel
##'
##' The 12 immune-checkpoint and cytotoxic-effector gene symbols used by
##' \code{plot_checkpoint_correlation} and \code{plot_gene_heatmap} when no
##' custom list is supplied.
##' @format A character vector of official gene symbols.
##' @name IMMUNE_CHECKPOINT_GENES
##' @export
IMMUNE_CHECKPOINT_GENES <- c("CD274", "PDCD1", "CD8A", "GZMA", "GZMB", "IFNG",
                             "CXCL9", "CXCL10", "PRF1", "TBX21", "CTLA4", "LAG3")

# ---------------------------------------------------------------------------
# internal helpers
# ---------------------------------------------------------------------------
.bc_colors <- function() {
  list(
    resp  = c("Non-responder" = "#2E9FDF", "Responder" = "#FF7777"),
    grp   = c("Low" = "#2E9FDF", "High" = "#FF7777"),
    pheno = c("desert" = "#2E9FDF", "excluded" = "#F5B400", "inflamed" = "#FF7777"),
    ic    = c("IC0" = "#2E9FDF", "IC1" = "#F5B400", "IC2+" = "#FF7777"),
    unknown = "#BFBFBF"
  )
}

.bc_no_grid <- ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                              panel.grid.minor = ggplot2::element_blank())

.bc_km_theme <- ggplot2::theme_classic() + ggplot2::theme(
  legend.title = ggplot2::element_text(size = 11),
  legend.text  = ggplot2::element_text(size = 10),
  plot.title   = ggplot2::element_text(hjust = 0.5, size = 13),
  axis.title   = ggplot2::element_text(size = 12),
  axis.text    = ggplot2::element_text(size = 11, color = "black")
)

.resp_factor <- function(cds) {
  r <- colData(cds)$binaryResponse
  factor(ifelse(r == "CR/PR", "Responder",
                ifelse(r == "SD/PD", "Non-responder", NA_character_)),
         levels = c("Non-responder", "Responder"))
}

.expr_group <- function(expr, gene) {
  x <- as.numeric(expr[gene, ])
  med <- median(x, na.rm = TRUE)
  factor(ifelse(x >= med, "High", "Low"), levels = c("Low", "High"))
}

.save_ggsurv <- function(g, outfile, width, height, dpi) {
  ext <- tolower(tools::file_ext(outfile))
  if (ext == "pdf") grDevices::pdf(outfile, width = width, height = height)
  else grDevices::png(outfile, width = width, height = height, units = "in", res = dpi)
  print(g)
  grDevices::dev.off()
  invisible(outfile)
}

# ---------------------------------------------------------------------------
#' Build a log2(TPM + 1) expression matrix from the IMvigor210 \code{cds}
#'
#' Computes transcripts-per-million (TPM) from raw counts and transcript length,
#' then applies the same edgeR TMM / CPM low-expression filter used by
#' \code{filterNvoom}. Returned values are \code{log2(TPM + 1)}.
#'
#' @param cds a \code{DESeqDataSet} (the bundled \code{cds} object).
#' @param genes optional character vector of gene symbols to subset to. Genes
#'   absent after filtering raise an error.
#' @param min_cpm minimum counts-per-million for the low-expression filter.
#' @param min_samples minimum number of samples that must reach \code{min_cpm};
#'   defaults to one tenth of the sample count (as in \code{filterNvoom}).
#' @param verbose logical; print filtering / cohort diagnostics.
#' @return A numeric matrix (genes x samples) of \code{log2(TPM + 1)} values with
#'   gene symbols as row names.
#' @export
get_tpm_expression <- function(cds, genes = NULL,
                               min_cpm = 0.25,
                               min_samples = NULL,
                               verbose = TRUE) {
  stopifnot(inherits(cds, "DESeqDataSet"))
  count_dat <- counts(cds)
  rd <- as.data.frame(rowData(cds))
  sym_col <- if (!is.null(rd$symbol)) "symbol" else "Symbol"
  eff <- rd[, c("entrez_id", "length", sym_col)]
  colnames(eff) <- c("entrez_id", "length", "symbol")
  eff <- eff[eff$symbol != "" & !is.na(eff$symbol), ]
  eff <- eff[order(eff$symbol, eff$length), ]   # same-name genes: keep shortest
  eff <- eff[!duplicated(eff$symbol), ]
  id <- intersect(rownames(eff), rownames(count_dat))
  eff <- eff[id, , drop = FALSE]
  count_dat <- count_dat[id, , drop = FALSE]
  rownames(count_dat) <- eff$symbol

  Counts2TPM <- function(counts, effLen) {
    rate  <- log(counts) - log(effLen)
    denom <- log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }
  tpm <- as.data.frame(apply(count_dat, 2, Counts2TPM, effLen = eff$length))

  dge <- DGEList(counts = count_dat)
  dge <- calcNormFactors(dge)
  cpm_mat <- cpm(dge, normalized.lib.sizes = TRUE)
  if (is.null(min_samples)) min_samples <- ncol(count_dat) / 10
  keep <- rowSums(cpm_mat >= min_cpm) >= min_samples
  if (verbose)
    message(sprintf("[get_tpm_expression] low-expression filter: %d/%d genes kept (minCpm=%.2f in >= %.0f samples)",
                    sum(keep), nrow(count_dat), min_cpm, min_samples))
  tpm <- tpm[keep, , drop = FALSE]
  expr <- log2(tpm + 1)

  if (!is.null(genes)) {
    miss <- setdiff(genes, rownames(expr))
    if (length(miss))
      stop("gene(s) not available after filtering: ", paste(miss, collapse = ", "))
    expr <- expr[genes, , drop = FALSE]
  }
  expr
}

# ---------------------------------------------------------------------------
#' Expression of a gene in responders vs. non-responders
#'
#' Violin + boxplot of a target gene's expression in the EVAL cohort, split by
#' binary response, with a Wilcoxon significance bracket.
#'
#' @param cds a \code{DESeqDataSet}.
#' @param gene character, gene symbol.
#' @param expr optional precomputed matrix from \code{get_tpm_expression}; if
#'   \code{NULL} it is computed internally.
#' @param outfile optional output path; the extension (\code{png}/\code{pdf})
#'   selects the format and a 2-panel combined figure is written.
#' @param width,height,dpi figure dimensions / resolution for \code{outfile}.
#' @param verbose logical; print the Wilcoxon p-value.
#' @return A list with \code{violin} and \code{boxplot} ggplot objects and the
#'   element \code{p.value} (Wilcoxon test).
#' @export
plot_expression_by_response <- function(cds, gene, expr = NULL,
                                        outfile = NULL,
                                        width = 3, height = 4, dpi = 300,
                                        verbose = TRUE) {
  if (is.null(expr)) expr <- get_tpm_expression(cds, verbose = verbose)
  stopifnot(gene %in% rownames(expr))
  resp <- .resp_factor(cds)
  df <- data.frame(Sample = droplevels(resp), exp = as.numeric(expr[gene, ]),
                   stringsAsFactors = FALSE)
  df <- df[!is.na(df$Sample), ]
  col <- .bc_colors()$resp
  cmp <- list(levels(df$Sample))

  box <- ggplot2::ggplot(df, ggplot2::aes(Sample, exp, fill = Sample)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = Sample), notch = FALSE,
                          position = ggplot2::position_dodge(width = 0.8),
                          outlier.alpha = 1, width = 0.4) +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = 0.2, add = 0.2)) +
    ggplot2::scale_fill_manual(values = col) +
    ggsignif::geom_signif(comparisons = cmp, step_increase = 0.1,
                          map_signif_level = TRUE, margin_top = 0.05,
                          test = "wilcox.test") +
    ggplot2::labs(y = paste0("Expression of ", gene, "\nlog2(TPM + 1)"),
                  title = "IMvigor210") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", colour = "black", linewidth = 0.25),
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.line = ggplot2::element_line(colour = "black", linewidth = 0.25),
      axis.title = ggplot2::element_text(size = 10, face = "plain", color = "black"),
      axis.text.x = ggplot2::element_text(face = "plain", colour = "black", vjust = 1),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(size = 12, face = "plain", color = "black"),
      axis.text = ggplot2::element_text(size = 12, color = "black"),
      legend.position = "none")

  violin <- ggplot2::ggplot(df, ggplot2::aes(x = Sample, y = exp)) +
    ggplot2::geom_violin(ggplot2::aes(fill = Sample), color = NA, alpha = 0.6,
                         width = 0.7, trim = TRUE, scale = "width") +
    ggplot2::geom_point(ggplot2::aes(color = Sample, fill = Sample),
                        show.legend = FALSE,
                        position = ggplot2::position_jitter(seed = 123456, width = 0.2),
                        shape = 21, size = 2) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = Sample), width = 0.5, size = 0.5,
                          alpha = 0.6, outlier.shape = NA) +
    ggsignif::geom_signif(comparisons = cmp, step_increase = 0.1,
                          map_signif_level = TRUE, margin_top = 0.2,
                          test = "wilcox.test") +
    ggplot2::scale_fill_manual(values = col) +
    ggplot2::scale_color_manual(values = col) +
    ggplot2::labs(y = paste0("Expression of ", gene, "\nlog2(TPM + 1)"),
                  title = "IMvigor210") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", colour = "black", linewidth = 0.25),
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.line = ggplot2::element_line(colour = "black", linewidth = 0.25),
      axis.title = ggplot2::element_text(size = 10, face = "plain", color = "black"),
      axis.text.x = ggplot2::element_text(face = "plain", colour = "black", vjust = 1),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(size = 12, face = "plain", color = "black"),
      axis.text = ggplot2::element_text(size = 12, color = "black"),
      legend.position = "none")

  p_w <- stats::wilcox.test(exp ~ Sample, data = df)$p.value
  if (verbose)
    message(sprintf("[plot_expression_by_response] %s: Responder vs Non-responder Wilcoxon p = %s (n=%d)",
                    gene, format.pval(p_w, digits = 3), nrow(df)))
  res <- list(violin = violin, boxplot = box, p.value = p_w)
  if (!is.null(outfile)) {
    ext <- tolower(tools::file_ext(outfile))
    combined <- violin + box + patchwork::plot_layout(ncol = 2)
    if (ext == "pdf") ggplot2::ggsave(outfile, combined, width = width * 2, height = height)
    else ggplot2::ggsave(outfile, combined, width = width * 2, height = height, dpi = dpi)
  }
  res
}

# ---------------------------------------------------------------------------
#' Correlation of a gene with immune-checkpoint genes
#'
#' Scatter "volcano" of the Pearson correlation coefficient (R) of a target gene
#' versus each gene in \code{checkpoints}, with significance highlighted.
#'
#' @param cds a \code{DESeqDataSet}.
#' @param gene character, gene symbol.
#' @param checkpoints character vector of gene symbols (default
#'   \code{IMMUNE_CHECKPOINT_GENES}).
#' @param expr optional precomputed matrix from \code{get_tpm_expression}.
#' @param outfile optional output path (\code{png}/\code{pdf}).
#' @param width,height,dpi figure dimensions / resolution for \code{outfile}.
#' @param verbose logical; print min correlation and number of significant genes.
#' @return The ggplot object; the per-gene result table is attached as the
#'   \code{table} attribute (\code{attr(p, "table")}).
#' @export
plot_checkpoint_correlation <- function(cds, gene,
                                        checkpoints = IMMUNE_CHECKPOINT_GENES,
                                        expr = NULL, outfile = NULL,
                                        width = 5.5, height = 5, dpi = 300,
                                        verbose = TRUE) {
  if (is.null(expr)) expr <- get_tpm_expression(cds, verbose = verbose)
  stopifnot(gene %in% rownames(expr))
  x <- as.numeric(expr[gene, ])
  rows <- lapply(checkpoints, function(g) {
    if (!g %in% rownames(expr)) return(NULL)
    y <- as.numeric(expr[g, ])
    ct <- stats::cor.test(x, y, method = "pearson")
    data.frame(checkpoint = g, R = ct$estimate, p = ct$p.value,
               neglogp = -log10(ct$p.value), absR = abs(ct$estimate))
  })
  res <- do.call(rbind, rows)
  col <- .bc_colors()
  LOW <- col$grp["Low"]; HIGH <- col$grp["High"]
  p <- ggplot2::ggplot(res, ggplot2::aes(x = R, y = neglogp)) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed",
                        color = "gray60", linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        color = "gray60", linewidth = 0.5) +
    ggplot2::geom_point(ggplot2::aes(size = absR, fill = R), shape = 21,
                        color = "black", stroke = 0.6, alpha = 0.85) +
    ggplot2::scale_size_continuous(range = c(3, 8), name = "|R|",
                                   guide = ggplot2::guide_legend(order = 2)) +
    ggplot2::scale_fill_gradient2(low = LOW, mid = "#F5F5F5", high = HIGH,
                                  midpoint = 0, name = "R",
                                  guide = ggplot2::guide_colorbar(order = 1)) +
    ggrepel::geom_text_repel(ggplot2::aes(label = checkpoint), size = 2.5,
                             color = "black", max.overlaps = 20, force = 3,
                             box.padding = 1.2, point.padding = 0.8,
                             min.segment.length = 0, segment.color = "gray50",
                             segment.size = 0.3) +
    ggplot2::labs(title = sprintf("%s correlation with immune checkpoint genes", gene),
                  x = "Correlation coefficient (R)",
                  y = expression(-log[10](P~value))) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 13, hjust = 0.5, face = "plain"),
      axis.title = ggplot2::element_text(size = 12, face = "plain", color = "black"),
      axis.text = ggplot2::element_text(size = 11, face = "plain", color = "black"),
      panel.background = ggplot2::element_rect(fill = "transparent", colour = "black"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      legend.position = "right",
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 9),
      plot.margin = ggplot2::margin(15, 15, 15, 15))
  if (verbose)
    message(sprintf("[plot_checkpoint_correlation] %s: min R=%.3f (vs %s), %d/%d checkpoints sig (p<0.05)",
                    gene, min(res$R), checkpoints[which.min(res$R)], sum(res$p < 0.05), nrow(res)))
  attr(p, "table") <- res
  if (!is.null(outfile)) {
    ext <- tolower(tools::file_ext(outfile))
    if (ext == "pdf") ggplot2::ggsave(outfile, p, width = width, height = height)
    else ggplot2::ggsave(outfile, p, width = width, height = height, dpi = dpi)
  }
  p
}

# ---------------------------------------------------------------------------
#' High/low expression groups and response proportion
#'
#' Percentage-stacked bar chart of responder / non-responder composition within
#' the high- and low-expression groups (median split) of a target gene, with a
#' chi-square test of independence.
#'
#' @param cds a \code{DESeqDataSet}.
#' @param gene character, gene symbol.
#' @param expr optional precomputed matrix from \code{get_tpm_expression}.
#' @param outfile optional output path (\code{png}/\code{pdf}).
#' @param width,height,dpi figure dimensions / resolution for \code{outfile}.
#' @param verbose logical; print the chi-square p-value.
#' @return The ggplot object; the chi-square p-value and the contingency table
#'   are attached as attributes \code{chisq_p} and \code{counts}.
#' @export
plot_response_proportion <- function(cds, gene, expr = NULL,
                                     outfile = NULL, width = 6, height = 5,
                                     dpi = 300, verbose = TRUE) {
  if (is.null(expr)) expr <- get_tpm_expression(cds, verbose = verbose)
  stopifnot(gene %in% rownames(expr))
  resp <- .resp_factor(cds)
  eg <- .expr_group(expr, gene)
  df <- data.frame(exprGroup = eg, response = droplevels(resp),
                   stringsAsFactors = FALSE)
  df <- df[!is.na(df$exprGroup) & !is.na(df$response), ]
  tab <- table(df$exprGroup, df$response)
  chi <- stats::chisq.test(tab)
  col <- .bc_colors()$resp
  prop <- as.data.frame(prop.table(tab, margin = 1))
  colnames(prop) <- c("exprGroup", "response", "pct")
  prop$pct <- prop$pct * 100
  prop$lab <- sprintf("%.1f%%", prop$pct)
  med <- median(as.numeric(expr[gene, ]), na.rm = TRUE)
  p <- ggplot2::ggplot(prop, ggplot2::aes(x = exprGroup, y = pct, fill = response)) +
    ggplot2::geom_col(position = "fill", width = 0.55, color = "black", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = lab),
                       position = ggplot2::position_fill(vjust = 0.5),
                       color = "white", size = 4.5, fontface = "bold") +
    ggplot2::scale_fill_manual(values = col) +
    ggplot2::labs(
      title = sprintf("%s high/low expression & response status", gene),
      subtitle = sprintf("chi-square p = %.2e (median split = %.3f)", chi$p.value, med),
      x = sprintf("%s expression group (median split)", gene),
      y = "Percentage of patients (%)") +
    ggplot2::scale_y_continuous(labels = function(x) sprintf("%d%%", x * 100)) +
    ggplot2::theme_bw() + .bc_no_grid +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 13),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10),
                   axis.text = ggplot2::element_text(size = 11),
                   axis.title = ggplot2::element_text(size = 12),
                   legend.title = ggplot2::element_blank())
  if (verbose)
    message(sprintf("[plot_response_proportion] %s: high/low x response chi-square p = %s (n=%d)",
                    gene, format.pval(chi$p.value, digits = 3), nrow(df)))
  attr(p, "chisq_p") <- chi$p.value
  attr(p, "counts") <- tab
  if (!is.null(outfile)) {
    ext <- tolower(tools::file_ext(outfile))
    if (ext == "pdf") ggplot2::ggsave(outfile, p, width = width, height = height)
    else ggplot2::ggsave(outfile, p, width = width, height = height, dpi = dpi)
  }
  p
}

# ---------------------------------------------------------------------------
#' Row-scaled expression heatmap of a gene set, ordered by response
#'
#' Row-scaled (z-score) heatmap of the supplied genes, with columns ordered by
#' response status and a response annotation bar.
#'
#' @param cds a \code{DESeqDataSet}.
#' @param genes character vector of gene symbols (target first is recommended).
#' @param expr optional precomputed matrix from \code{get_tpm_expression}.
#' @param response_only logical; if \code{TRUE} keep only response-known samples
#'   (EVAL cohort), otherwise include all samples with \code{Unknown} responses
#'   shown in grey (ALL cohort).
#' @param outfile optional output path (\code{png}/\code{pdf}). When omitted the
#'   heatmap is drawn on the current device and the pheatmap object is returned.
#' @param width figure width in inches; height auto-scales with gene count.
#' @param dpi resolution for a \code{png} \code{outfile}.
#' @param verbose logical; print canvas dimensions.
#' @return Invisibly the output path (when \code{outfile} is given) or the
#'   \pkg{pheatmap} object (otherwise).
#' @export
plot_gene_heatmap <- function(cds, genes, expr = NULL, response_only = FALSE,
                              outfile = NULL, width = 8, dpi = 300,
                              verbose = TRUE) {
  if (is.null(expr)) expr <- get_tpm_expression(cds, verbose = verbose)
  miss <- setdiff(genes, rownames(expr))
  if (length(miss)) stop("gene(s) not in expression matrix: ", paste(miss, collapse = ", "))
  mat <- expr[genes, , drop = FALSE]
  resp <- .resp_factor(cds)
  col <- .bc_colors()
  CELL_H <- 14
  heat_colors <- grDevices::colorRampPalette(c(col$grp["Low"], "#F5F5F5", col$grp["High"]))(100)
  if (response_only) {
    keep <- which(!is.na(resp))
    ord <- keep[order(resp[keep])]
    mat <- mat[, ord, drop = FALSE]
    anno <- data.frame(Response = factor(as.character(resp[ord]),
                                         levels = c("Non-responder", "Responder")))
    lv <- c("Non-responder", "Responder")
  } else {
    ord <- order(resp)
    mat <- mat[, ord, drop = FALSE]
    anno <- data.frame(Response = factor(ifelse(is.na(resp[ord]), "Unknown",
                                                as.character(resp[ord])),
                                         levels = c("Non-responder", "Responder", "Unknown")))
    lv <- c("Non-responder", "Responder", "Unknown")
  }
  rownames(anno) <- colnames(cds)[ord]
  acl <- list(Response = c("Non-responder" = col$resp[["Non-responder"]],
                           "Responder" = col$resp[["Responder"]],
                           "Unknown" = col$unknown)[lv])
  heat_height <- nrow(mat) * CELL_H / 72 + 1.6
  args <- list(mat = mat, scale = "row", color = heat_colors,
               annotation_col = anno, annotation_colors = acl,
               cluster_cols = FALSE, cluster_rows = FALSE,
               show_rownames = TRUE, show_colnames = FALSE,
               cellheight = CELL_H, fontsize_row = 10, fontsize_col = 8,
               main = sprintf("%s (row-scaled, %s n=%d)",
                              paste(genes, collapse = "/"),
                              ifelse(response_only, "response-known", "all"), ncol(mat)),
               width = width, height = heat_height)

  if (is.null(outfile)) {
    return(do.call(pheatmap::pheatmap, args))
  }
  ext <- tolower(tools::file_ext(outfile))
  if (ext == "pdf") {
    tmp <- tempfile(fileext = ".pdf")
    do.call(pheatmap::pheatmap, c(args, list(filename = tmp)))
    if (file.exists(outfile)) {
      ok <- try(file.remove(outfile), silent = TRUE)
      if (inherits(ok, "try-error") || file.exists(outfile)) {
        alt <- sub("\\.pdf$", "_plot.pdf", outfile)
        file.rename(tmp, alt); outfile <- alt
      } else file.rename(tmp, outfile)
    } else file.rename(tmp, outfile)
  } else {
    do.call(pheatmap::pheatmap, c(args, list(filename = outfile)))
  }
  if (verbose)
    message(sprintf("[plot_gene_heatmap] %d genes, canvas height = %.2f in, response_only=%s",
                    nrow(mat), heat_height, response_only))
  invisible(outfile)
}

# ---------------------------------------------------------------------------
#' ROC / AUC for response prediction by a single gene
#'
#' Receiver-operating-characteristic curve of a target gene's expression for
#' predicting responder status, with the AUC and 95\% CI.
#'
#' @param cds a \code{DESeqDataSet}.
#' @param gene character, gene symbol.
#' @param expr optional precomputed matrix from \code{get_tpm_expression}.
#' @param outfile optional output path (\code{png}/\code{pdf}).
#' @param width,height,dpi figure dimensions / resolution for \code{outfile}.
#' @param verbose logical; print AUC and 95\% CI.
#' @return The ggplot object; the AUC value and its CI are attached as
#'   attributes \code{auc} and \code{ci}.
#' @export
plot_roc_response <- function(cds, gene, expr = NULL, outfile = NULL,
                              width = 5, height = 5, dpi = 300, verbose = TRUE) {
  if (is.null(expr)) expr <- get_tpm_expression(cds, verbose = verbose)
  stopifnot(gene %in% rownames(expr))
  resp <- .resp_factor(cds)
  df <- data.frame(gene = as.numeric(expr[gene, ]), response = droplevels(resp),
                   stringsAsFactors = FALSE)
  df <- df[!is.na(df$response), ]
  resp01 <- ifelse(df$response == "Responder", 1, 0)
  rocobj <- pROC::roc(resp01 ~ df$gene, data = df)
  auc_val <- as.numeric(pROC::auc(rocobj))
  ci <- as.numeric(pROC::ci.auc(rocobj))
  col <- .bc_colors()
  HIGH <- col$grp["High"]
  roc_df <- data.frame(fpr = 1 - rocobj$specificities, tpr = rocobj$sensitivities)
  roc_df <- roc_df[order(roc_df$fpr, roc_df$tpr), ]
  roc_df <- rbind(c(0, 0), roc_df, c(1, 1))
  p <- ggplot2::ggplot(roc_df, ggplot2::aes(x = fpr, y = tpr)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         color = "gray60", linewidth = 0.5) +
    ggplot2::geom_line(color = HIGH, linewidth = 1) +
    ggplot2::annotate("text", x = 0.55, y = 0.12,
                      label = sprintf("AUC = %.3f\n(95%% CI %.3f-%.3f)", auc_val, ci[1], ci[3]),
                      size = 4, hjust = 0) +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, 0.25)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, 0.25)) +
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::labs(title = sprintf("%s predicts response (ROC)", gene),
                  x = "1 - Specificity", y = "Sensitivity") +
    ggplot2::theme_bw() + .bc_no_grid +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 13),
                   axis.title = ggplot2::element_text(size = 12),
                   axis.text  = ggplot2::element_text(size = 11, color = "black"))
  if (verbose)
    message(sprintf("[plot_roc_response] %s: AUC = %.3f (95%% CI %.3f-%.3f), n=%d",
                    gene, auc_val, ci[1], ci[3], nrow(df)))
  attr(p, "auc") <- auc_val
  attr(p, "ci") <- ci
  if (!is.null(outfile)) {
    ext <- tolower(tools::file_ext(outfile))
    if (ext == "pdf") ggplot2::ggsave(outfile, p, width = width, height = height)
    else ggplot2::ggsave(outfile, p, width = width, height = height, dpi = dpi)
  }
  p
}

# ---------------------------------------------------------------------------
#' Kaplan-Meier survival by expression level or by response
#'
#' Overall-survival Kaplan-Meier curve, stratified either by high/low expression
#' (median split, ALL cohort) or by binary response (EVAL cohort).
#'
#' @param cds a \code{DESeqDataSet}.
#' @param gene character, gene symbol.
#' @param by character; \code{"exprGroup"} (default) stratifies by high/low
#'   expression, \code{"response"} stratifies by responder status.
#' @param expr optional precomputed matrix from \code{get_tpm_expression}.
#' @param outfile optional output path (\code{png}/\code{pdf}).
#' @param width,height,dpi figure dimensions / resolution for \code{outfile}.
#' @param verbose logical; print the log-rank p-value.
#' @return A \pkg{survminer} \code{ggsurvplot} object; the log-rank p-value is
#'   attached as the \code{pvalue} attribute.
#' @export
plot_survival <- function(cds, gene, by = c("exprGroup", "response"),
                          expr = NULL, outfile = NULL, width = 5, height = 5,
                          dpi = 300, verbose = TRUE) {
  by <- match.arg(by)
  if (is.null(expr)) expr <- get_tpm_expression(cds, verbose = verbose)
  stopifnot(gene %in% rownames(expr))
  col <- .bc_colors()
  os <- colData(cds)$os
  cens <- colData(cds)$censOS
  if (by == "exprGroup") {
    grp <- .expr_group(expr, gene)
    df <- data.frame(exprGroup = grp, os = os, censOS = cens, stringsAsFactors = FALSE)
    df <- df[!is.na(df$exprGroup) & !is.na(df$censOS), ]
    fit <- survival::survfit(survival::Surv(os, censOS) ~ exprGroup, data = df)
    pv <- survminer::surv_pvalue(fit, data = df)$pval[1]
    g <- survminer::ggsurvplot(fit, data = df, pval = TRUE, pval.coord = c(0, 0.05),
                               palette = c(col$grp["Low"], col$grp["High"]),
                               xlab = "Time (months)", ylab = "Overall survival probability",
                               legend.title = gene, legend.labs = c("Low", "High"),
                               risk.table = FALSE, ggtheme = .bc_km_theme,
                               title = sprintf("%s high/low & overall survival", gene))
    if (verbose) message(sprintf("[plot_survival] %s by expression: log-rank p = %s (n=%d)",
                                 gene, format.pval(pv, digits = 3), nrow(df)))
    attr(g, "pvalue") <- pv
  } else {
    resp <- .resp_factor(cds)
    df <- data.frame(response = droplevels(resp), os = os, censOS = cens, stringsAsFactors = FALSE)
    df <- df[!is.na(df$response) & !is.na(df$censOS), ]
    fit <- survival::survfit(survival::Surv(os, censOS) ~ response, data = df)
    pv <- survminer::surv_pvalue(fit, data = df)$pval[1]
    g <- survminer::ggsurvplot(fit, data = df, pval = TRUE, pval.coord = c(0, 0.05),
                               palette = c(col$resp["Non-responder"], col$resp["Responder"]),
                               xlab = "Time (months)", ylab = "Overall survival probability",
                               legend.title = gene, legend.labs = c("Non-responder", "Responder"),
                               risk.table = FALSE, ggtheme = .bc_km_theme,
                               title = sprintf("%s: survival by response", gene))
    if (verbose) message(sprintf("[plot_survival] %s by response: log-rank p = %s (n=%d)",
                                 gene, format.pval(pv, digits = 3), nrow(df)))
    attr(g, "pvalue") <- pv
  }
  if (!is.null(outfile)) .save_ggsurv(g, outfile, width, height, dpi)
  g
}

# ---------------------------------------------------------------------------
# Multi-group violin helper (internal styling)
# ---------------------------------------------------------------------------
.violin_multi <- function(d, gene, colors, title) {
  d <- d[!is.na(d$group) & !is.na(d$value), ]
  d$group <- droplevels(d$group)
  lv <- levels(d$group)
  cmp <- combn(lv, 2, simplify = FALSE)
  kw <- stats::kruskal.test(value ~ group, data = d)$p.value
  ggplot2::ggplot(d, ggplot2::aes(x = group, y = value, fill = group)) +
    ggplot2::geom_violin(color = NA, alpha = 0.6, width = 0.7, trim = TRUE, scale = "width") +
    ggplot2::geom_point(ggplot2::aes(fill = group), show.legend = FALSE,
                        position = ggplot2::position_jitter(seed = 123456, width = 0.2),
                        shape = 21, size = 2, color = "black") +
    ggplot2::geom_boxplot(width = 0.5, size = 0.5, alpha = 0.6, outlier.shape = NA) +
    ggsignif::geom_signif(comparisons = cmp, step_increase = 0.11, map_signif_level = TRUE,
                          margin_top = 0.08, tip_length = 0.01, test = "wilcox.test",
                          textsize = 3.5) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.16))) +
    ggplot2::labs(y = paste0("Expression of ", gene, "\nlog2(TPM + 1)"),
                  x = NULL, title = title,
                  subtitle = sprintf("Kruskal-Wallis p = %s", format.pval(kw, digits = 3))) +
    ggplot2::theme_classic() +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", colour = "black", linewidth = 0.25),
                   axis.line = ggplot2::element_line(colour = "black", linewidth = 0.25),
                   axis.title = ggplot2::element_text(size = 10, face = "plain", color = "black"),
                   axis.text = ggplot2::element_text(size = 11, color = "black"),
                   axis.title.y = ggplot2::element_text(size = 12, face = "plain", color = "black"),
                   plot.title = ggplot2::element_text(hjust = 0.5, size = 12),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9, color = "gray30"),
                   legend.position = "none")
}

# ---------------------------------------------------------------------------
#' Immune-microenvironment context of a gene's expression
#'
#' Violin plots of a target gene across immune phenotypes (desert / excluded /
#' inflamed) and PD-L1 immune-cell (IC) levels (IC0 / IC1 / IC2+), in the IMM
#' cohort, with pairwise Wilcoxon brackets and a Kruskal-Wallis p-value.
#'
#' @param cds a \code{DESeqDataSet}.
#' @param gene character, gene symbol.
#' @param expr optional precomputed matrix from \code{get_tpm_expression}.
#' @param outfile optional output path (\code{png}/\code{pdf}); a 2-panel
#'   (phenotype over IC level) figure is written.
#' @param width,height,dpi figure dimensions / resolution for \code{outfile}.
#' @param verbose logical; print the Kruskal-Wallis p-value.
#' @return A list with \code{phenotype} and \code{ic} ggplot objects and the
#'   element \code{kruskal_p}.
#' @export
plot_immune_context <- function(cds, gene, expr = NULL, outfile = NULL,
                                width = 4.2, height = 5.5, dpi = 300, verbose = TRUE) {
  if (is.null(expr)) expr <- get_tpm_expression(cds, verbose = verbose)
  stopifnot(gene %in% rownames(expr))
  cd <- as.data.frame(colData(cds))
  keep <- !is.na(cd$IC.Level) & !is.na(cd$binaryResponse) & !is.na(cd$Immune.phenotype)
  cd <- cd[keep, ]
  val <- as.numeric(expr[gene, colnames(cds)[keep]])
  col <- .bc_colors()
  d_pheno <- data.frame(group = factor(cd$Immune.phenotype, levels = names(col$pheno)),
                        value = val)
  p_pheno <- .violin_multi(d_pheno, gene, col$pheno,
                           sprintf("%s by immune phenotype", gene))
  d_ic <- data.frame(group = factor(cd$IC.Level, levels = names(col$ic)),
                     value = val)
  p_ic <- .violin_multi(d_ic, gene, col$ic,
                        sprintf("%s by immune cell (IC) level", gene))
  kw_p <- stats::kruskal.test(val ~ factor(cd$Immune.phenotype))$p.value
  if (verbose)
    message(sprintf("[plot_immune_context] %s by immune phenotype Kruskal p = %s (n=%d)",
                    gene, format.pval(kw_p, digits = 3), length(val)))
  res <- list(phenotype = p_pheno, ic = p_ic, kruskal_p = kw_p)
  if (!is.null(outfile)) {
    ext <- tolower(tools::file_ext(outfile))
    combined <- p_pheno + p_ic + patchwork::plot_layout(ncol = 1)
    if (ext == "pdf") ggplot2::ggsave(outfile, combined, width = width, height = height)
    else ggplot2::ggsave(outfile, combined, width = width, height = height, dpi = dpi)
  }
  res
}

# ---------------------------------------------------------------------------
#' Run the full biomarker analysis for one gene and write all figures
#'
#' Convenience wrapper that executes every helper above on a single gene and
#' writes the figures and summary tables into a folder structure mirroring the
#' original standalone analysis script.
#'
#' @param cds a \code{DESeqDataSet}.
#' @param gene character gene symbol.
#' @param outdir output directory (created if missing).
#' @param checkpoints character vector (default \code{IMMUNE_CHECKPOINT_GENES}).
#' @param verbose logical.
#' @return A list with \code{plots} (named list of plot objects) and
#'   \code{summary} (data.frame of key metrics, also written to
#'   \code{<outdir>/tables/<gene>_summary_metrics.csv}).
#' @export
run_biomarker_report <- function(cds, gene,
                                 outdir = file.path("Results", gene),
                                 checkpoints = IMMUNE_CHECKPOINT_GENES,
                                 verbose = TRUE) {
  expr <- get_tpm_expression(cds, verbose = verbose)
  stopifnot(gene %in% rownames(expr))
  dirs <- c("expr_diff", "correlation", "group_proportion", "heatmap",
            "survival", "roc", "immune_context", "tables")
  for (d in dirs) dir.create(file.path(outdir, d), recursive = TRUE, showWarnings = FALSE)

  resp <- .resp_factor(cds)
  N_ALL <- ncol(cds)
  R_ALL <- sum(resp == "Responder", na.rm = TRUE)
  N_EVAL <- sum(!is.na(resp))
  R_EVAL <- sum(resp == "Responder", na.rm = TRUE)
  cd <- as.data.frame(colData(cds))
  imm_keep <- !is.na(cd$IC.Level) & !is.na(cd$binaryResponse) & !is.na(cd$Immune.phenotype)
  N_IMM <- sum(imm_keep)
  R_IMM <- sum(resp[imm_keep] == "Responder", na.rm = TRUE)

  plots <- list()

  p1 <- plot_expression_by_response(cds, gene, expr = expr, verbose = verbose)
  ggplot2::ggsave(file.path(outdir, "expr_diff", sprintf("%s_violin_resp%d.pdf", gene, N_EVAL)),
                  p1$violin + p1$boxplot + patchwork::plot_layout(ncol = 2), width = 6, height = 4)
  plots$expression_by_response <- p1

  p2 <- plot_checkpoint_correlation(cds, gene, checkpoints = checkpoints, expr = expr, verbose = verbose)
  ggplot2::ggsave(file.path(outdir, "correlation", sprintf("%s_corr.pdf", gene)), p2, width = 5.5, height = 5)
  utils::write.csv(attr(p2, "table"), file.path(outdir, "tables", sprintf("%s_correlation_table.csv", gene)), row.names = FALSE)
  plots$checkpoint_correlation <- p2

  p3 <- plot_response_proportion(cds, gene, expr = expr, verbose = verbose)
  ggplot2::ggsave(file.path(outdir, "group_proportion", sprintf("%s_resp_proportion.pdf", gene)), p3, width = 6, height = 5)
  plots$response_proportion <- p3

  genes_h <- c(gene, checkpoints)
  plot_gene_heatmap(cds, genes_h, expr = expr, response_only = FALSE,
                    outfile = file.path(outdir, "heatmap", sprintf("%s_heatmap_all%d.pdf", gene, N_ALL)), verbose = verbose)
  plot_gene_heatmap(cds, genes_h, expr = expr, response_only = TRUE,
                    outfile = file.path(outdir, "heatmap", sprintf("%s_heatmap_resp%d.pdf", gene, N_EVAL)), verbose = verbose)

  p5a <- plot_survival(cds, gene, by = "exprGroup", expr = expr, verbose = verbose)
  .save_ggsurv(p5a, file.path(outdir, "survival", sprintf("%s_KM_highlow_all%d.pdf", gene, N_ALL)), 5, 5, 300)
  p5b <- plot_survival(cds, gene, by = "response", expr = expr, verbose = verbose)
  .save_ggsurv(p5b, file.path(outdir, "survival", sprintf("%s_KM_response_resp%d.pdf", gene, N_EVAL)), 5, 5, 300)
  plots$survival_exprGroup <- p5a
  plots$survival_response <- p5b

  p6 <- plot_roc_response(cds, gene, expr = expr, verbose = verbose)
  ggplot2::ggsave(file.path(outdir, "roc", sprintf("%s_ROC_resp%d.pdf", gene, N_EVAL)), p6, width = 5, height = 5)
  plots$roc <- p6

  p7 <- plot_immune_context(cds, gene, expr = expr, verbose = verbose)
  ggplot2::ggsave(file.path(outdir, "immune_context", sprintf("%s_immune_context_resp%d.pdf", gene, N_IMM)),
                  p7$phenotype + p7$ic + patchwork::plot_layout(ncol = 1), width = 4.2, height = 5.5)
  plots$immune_context <- p7

  summary <- data.frame(
    metric = c("cohort_ALL_n", "cohort_ALL_responder", "cohort_EVAL_n", "cohort_EVAL_responder",
               "cohort_IMMUNE_n", "cohort_IMMUNE_responder",
               "wilcox_Responder_vs_NonResponder_p",
               "min_corr_R_with_checkpoints", "n_checkpoints_sig_p05",
               "chi_sq_highlow_vs_response_p", "AUC_predict_response",
               "AUC_95CI_low", "AUC_95CI_high",
               "KM_logrank_highlow_p", "KM_logrank_response_p",
               "Kruskal_p_by_immune_phenotype"),
    value = c(N_ALL, R_ALL, N_EVAL, R_EVAL, N_IMM, R_IMM,
              p1$p.value,
              min(attr(p2, "table")$R), sum(attr(p2, "table")$p < 0.05),
              attr(p3, "chisq_p"), attr(p6, "auc"),
              attr(p6, "ci")[1], attr(p6, "ci")[3],
              attr(p5a, "pvalue"), attr(p5b, "pvalue"), p7$kruskal_p),
    stringsAsFactors = FALSE)
  utils::write.csv(summary, file.path(outdir, "tables", sprintf("%s_summary_metrics.csv", gene)), row.names = FALSE)
  if (verbose) message(sprintf("[run_biomarker_report] wrote all figures + tables to %s", outdir))
  list(plots = plots, summary = summary)
}
