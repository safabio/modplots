#' Custom pheatmap of DESeq2 vst data
#'
#' Preset heatmap for DESeq2 data. Plots z scored (scaled to mean = 0 and sd = 1) vst expression of any gene list.
#'
#' @param vsd vsd object (created as following: \code{vsd <- vst(dds)})
#' @param IDs vector of gene ids
#' @param ... parameters passed to pheatmap
#'
#' @importFrom pheatmap pheatmap
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom dplyr left_join
#' @importFrom grDevices colorRampPalette
#'
#' @export
#'
#'

mPheatmapDESeq2 <- function(vsd,
                            IDs,
                            ...
                            ) {

    gnames <- modplots::gnames

    heat_col <- colorRampPalette(colors = c("darkblue","dodgerblue4", "white", "red", "darkred"))

    toplot <- assay(vsd)[IDs, ]
    toplot <- t(scale(t(toplot)))
    toplot <- as.data.frame(toplot)
    toplot <- rownames_to_column(toplot, "Gene.stable.ID")
    toplot <- left_join(toplot, gnames, by = "Gene.stable.ID")
    toplot <- column_to_rownames(toplot, "Gene.name")

    toplot$Gene.stable.ID <- NULL

    lim <- max(floor(min(toplot)*100), floor(max(toplot)*100))

    FDQ_AM_lig_htmp <- pheatmap(
      toplot,
      color = heat_col(2*lim)[lim+(floor(min(toplot)*100):floor(max(toplot)*100))],
      fontsize = 7,
      border_color = NA,
      ...
    )

}


