#' Volcano plot visualization
#'
#' Intuitive way of visualizing log fold changes and adjusted pvalues
#'
#' @param table DE table created by either DESeq2 or Seurat.
#' For DESeq2 the table need the columns ID, log2FoldChange, and padj.
#' The Seurat table needs the colums gene, avg_log2_FC, and p_val_adj.
#' @param type character indicating the method that created the DE table. Can be either seurat, deseq2 or none. If "none" is
#' specified the function looks for "gene" or "id", "adj", and "2F" in the tables colnames.
#' @param l2fc log2 fold change threshold
#' @param p.adj Adjusted p-value threshold
#' @param return Boolean whether to return a plot (FALSE, default) or a ggplot objects (TRUE)
#' @param plotly Boolean whether to plot the volcano plot with plotly to create an interactive plot.
#' @param gnames Optional, data frame of two columns, Gene.stable.ID and Gene.name. If none is specified, internal gg6 chicken gnames table is used.
#'
#' @details Adjusted pvals of 0 are capped to 1e-350 so they are not driven to infinity by -log10 transformation.
#'
#' @importFrom dplyr mutate select left_join case_when
#' @importFrom ggplot2 geom_point geom_hline geom_vline scale_fill_manual scale_size_manual scale_alpha_manual ylab
#' theme_bw aes
#' @importFrom plotly ggplotly
#'
#' @return A ggplot plot if \code{return = FALSE}, or else a ggplot object.
#'
#' @export
#'

mVolcanoplot <- function(table,
                         type = "seurat",
                         l2fc = log2(2),
                         p.adj = 0.05,
                         return = FALSE,
                         plotly = FALSE,
                         gnames = NULL) {
  if (!type %in% c("seurat", "deseq2", "none")) {
    stop("type should be either seurat, deseq2, or none")
  }

  if (is.null(gnames)) {
    gnames <- modplots::gnames
  } else if (!all(colnames(gnames) %in% c("Gene.name", "Gene.stable.ID"))) {
    stop("gnames should contain two rows and two rows only, called: Gene.stable.ID and Gene.name")
  }

  if (type == "seurat") {
    table <- select(table, gene, avg_log2FC, p_val_adj)
    colnames(table) <- c("Gene.stable.ID", "log2FC", "padj")
    table <- left_join(table, gnames, by = "Gene.stable.ID")
  } else if (type == "deseq2") {
    table <- select(table, ID, log2FoldChange, padj)
    colnames(table) <- c("Gene.stable.ID", "log2FC", "padj")
    table <- left_join(table, gnames, by = "Gene.stable.ID")
  } else if (type == "none") {
    id <- grep("gene|id", colnames(table),ignore.case = TRUE, value = TRUE)
    pa <- grep("adj", colnames(table),ignore.case = TRUE, value = TRUE)
    fc <- grep("2F", colnames(table),ignore.case = TRUE, value = TRUE)
    table <- select(table, c(id, fc, pa))
    colnames(table) <- c("Gene.stable.ID", "log2FC", "padj")
    table <- left_join(table, gnames, by = "Gene.stable.ID")
  }

  # Add colour, size and alpha (transparency) to volcano plot
  cols <- c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey")
  sizes <- c("up" = 2, "down" = 2, "ns" = 1)
  alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)

  toplot <-  mutate(table, gene_type = case_when(log2FC >= l2fc & padj <= p.adj ~ "up",
                                 log2FC <= -l2fc & padj <= p.adj ~ "down",
                                 TRUE ~ "ns"))

  toplot$padj <- -log10(toplot$padj)
  toplot$padj[is.infinite(toplot$padj)] <- 350


  volplot <-   ggplot(data = toplot,
                     aes(x = log2FC,
                         y = padj,
                         label = Gene.name,
                         fill = gene_type,
                         size = gene_type,
                         alpha = gene_type)) +
    geom_point(shape = 21) +
    geom_hline(yintercept = -log10(p.adj), linetype = "dashed") +
    geom_vline(xintercept = c(-l2fc,l2fc), linetype = "dashed") +
    scale_fill_manual(values = cols) + # Modify point colour
    scale_size_manual(values = sizes) + # Modify point size
    scale_alpha_manual(values = alphas) + # Modify point transparency
    ylab("-log10(padj)") +
    theme_bw()

  if(return) {
    return(volplot)
  } else {
    if (plotly) {
      ggplotly(volplot)
    } else {
      volplot
    }
  }

}

