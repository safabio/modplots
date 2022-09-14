#' Visualize 'features' on a dimensional reduction plot
#'
#' Colors single cells on a dimensional reduction plot according to a 'feature'
#' (i.e. gene expression, PC scores, number of genes detected, etc.)
#'
#' @param my.se Seurat object
#' @param my.reduct name of reduction to take from the embedding slot. e.g.: pca, tsne, umap ...
#' @param my.slot slot to take the data from. e.g.: counts, data, scale.data
#' @param my.features Vector of features to plot. Features must be present in modplots::gnames$Gene.name!
#' @param colors Vector of colors, used in scale_colour_gradientn()
#' @param size pt size of geom_point
#' @param alpha alpha parameter for geom_point
#' @param order Boolean determining whether to plot cells in order of expression. Can be useful if
#' cells expressing given feature are getting buried.
#' @param return Boolean, wether to return a plot (FALSE, default) or a list of ggplot objects (TRUE)
#'
#' @return A grid.arrange plot if
#' \code{return = FALSE}; otherwise, a list of ggplot objects
#'
#'
#' @import grDevices
#' @import Seurat
#' @import rlang
#' @import gridExtra
#' @import RColorBrewer
#' @import ggplot2
#'
#' @export
#'
#'
#'
#'
#'

mFeaturePlot <- function(my.se,
                         my.reduct = "tsne",
                         my.slot = "data",
                         my.features = NULL,
                         colors = c("gray90","gray90","gray80","yellow", "orange", "red", "darkred", "darkred"),
                         size=2,
                         alpha = 0.4,
                         order = TRUE,
                         return = FALSE) {

  if (!any(class(my.se) == "Seurat")) {
    stop("my.se must be an object of class Seurat")
  }

  emb <- data.frame(Embeddings(my.se, my.reduc))
  colnames(emb) <- c("reduc_1", "reduc_2")

  gnames <- modplots::gnames

  my.features <- toupper(my.features)

  if(!all(my.features %in% gnames$Gene.name)) {
    stop("Some or all features not found in gnames$Gene.name!")
  }

  GOI <- gnames[gnames$Gene.name %in% my.features,]
  GOI <- GOI[match(my.features, GOI$Gene.name),]

  if(!all(GOI$Gene.stable.ID %in% rownames(my.se))) {
    toremove <- GOI$Gene.name[!GOI$Gene.stable.ID %in% rownames(my.se)]

    GOI <- GOI[GOI$Gene.stable.ID %in% rownames(my.se),]

    cat(c("Not all features are present in your object! Removing: ", toremove, "\n"))
  }


  expr <- t(as.matrix(GetAssayData(my.se, slot =  my.slot)[GOI$Gene.stable.ID,]))

  if(!identical(rownames(expr), rownames(emb))) {
    stop("Embeddings and Assay do not contain the same cells!")
  }

  my.plots = list()

  for (i in 1:(nrow(GOI))) {
    ID <- GOI[i, "Gene.stable.ID"]
    name <- GOI[i, "Gene.name"]

    tmp <- cbind(emb, expr[,i])
    colnames(tmp) <- c("reduc_1","reduc_2","gene")

    if(order) {
      tmp <- dplyr::arrange(tmp, gene)
    }

    my.plots[[i]] <- ggplot(tmp, aes(x = reduc_1, y = reduc_2, color = gene)) +
      geom_point(size=size, alpha = 0.4, pch = 19) +
      scale_colour_gradientn(colours = colors) +
      theme_classic() + labs(colour=name) + ggtitle(name) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  }

  if (return) {
    return(my.plots)
  } else {

  grid.arrange(grobs = my.plots, ncol = ceiling(sqrt(nrow(GOI))))

  }
}
