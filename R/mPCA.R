#' mPCA
#'
#' Function to plot PCA biplot from a DESeq2 vst (variance stabilising transformation).
#'
#' @param vsd Variance stabilising transformed dds (vsd <- vst(dds)).
#' @param ntop integer, number of variable genes to use for PCA.
#' @param PCs character vector of length two, containing the two PCs to plot (default to: c("PC1", "PC2")).
#' @param genes vector of gene IDs to use instead of internally calculated variable genes.
#' @param meta Meta data frame for ggplot aestetics, if not specified grouping will be taken from \code{vsd`@`colData}.
#' If you specify a meta data frame, you have to also specify group, which should be a column of the meta data frame.
#' @param group Column of vsd`@`colData used as aesthetic for grouping the PCA. Meta data column if meta data is specified.
#' @param label Column of vsd`@`colData used as aesthetic for labeling the PCA. Meta data column if meta data is specified.
#' @param shape Column of vsd`@`colData used as aesthetic for point shapes in the PCA. Meta data column if meta data is specified.
#' @param pch vector of point characters if shape is specified, same length as factor levels in shape.
#' @param colors color vector, should be same length as factor levels in group
#' @param return Boolean whether to return a plot (Default = FALSE) or a ggplot object list (TRUE)
#'
#' @importFrom MatrixGenerics rowVars
#' @importFrom SummarizedExperiment assay
#' @importFrom stats prcomp
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join
#' @importFrom ggplot2 ggplot geom_point scale_shape_manual xlab ylab scale_color_manual
#' @importFrom BiocGenerics as.data.frame
#' @importFrom ggrepel geom_label_repel
#'
#' @export
#'


mPCA <- function(vsd,
                 ntop,
                 PCs = c("PC1", "PC2"),
                 genes = NULL,
                 meta = NULL,
                 group = NULL,
                 label = NULL,
                 shape = NULL,
                 pch = NULL,
                 colors = NULL,
                 return = FALSE) {

  if (!class(vsd)[[1]] == "DESeqTransform") {
    stop("vsd is not of class DESeqTransform!")
  }

  if (is.null(group)) {
    stop("You must specify a grouping variable!")
  }

  if (!is.null(meta)) {

    if (!is.data.frame(meta)) {
      stop("meta should be a data frame!")
    }

    if (!group %in% colnames(meta)) {
      stop("group should be a column of meta!")
    }

    if (!is.null(label)) {
      if (!label %in% colnames(meta)) {
      stop("label should be a column of meta!")
      }
    }

    if (!is.null(shape)) {
      if (!shape %in% colnames(meta)) {
        stop("shape should be a column of meta!")
      }
    }

    # if (!is.null(shape) & !is.null(pch) & !identical(length(pch), length(levels(meta[, shape])))) {
    #   stop("pch should be the same lenght as levels in the meta[, shape] column!")
    # }
  }


  ntop <- as.integer(ntop)

  if (!is.integer(ntop)) {
    stop("ntop must be of type integer!")
  }


  if (is.null(genes)) {
    rv <- rowVars(assay(vsd))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  } else {
    select <- genes
  }

  pca <- prcomp(t(assay(vsd)[select, ]))

  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

  xvar <- paste0(PCs[1], " : ",round(percentVar[1] * 100),"% variance")
  yvar <- paste0(PCs[2], " : ",round(percentVar[2] * 100),"% variance")

  PC <- as.data.frame(pca$x)
  PC <- rownames_to_column(PC, "sample")

  if (!is.null(meta)) {
    PC <- left_join(PC, meta, by = "sample")
  } else {
    coldat <- as.data.frame(vsd@colData)
    coldat <- rownames_to_column(coldat, "sample")
    PC <- left_join(PC, coldat, by = "sample")
  }


  PC_plot <- ggplot(data = PC, aes_string(x="PC1", y="PC2", color = group)) +
    geom_point(size=3) +
    xlab(xvar) +
    ylab(yvar)

  if (!is.null(shape)) {
    PC_plot <- PC_plot + scale_shape_manual(aes_string(shape = shape), values = pch)
  }

  if (!is.null(label)) {
    PC_plot <- PC_plot + geom_label_repel(aes_string(label = label))
  }

  if (!is.null(colors)) {
    PC_plot <- PC_plot + scale_color_manual(values = colors)
  }

 if (return == TRUE) {
   return(PC_plot)
 } else {
   PC_plot
 }

}

