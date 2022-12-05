#' loadingPlot
#'
#' Takes the rotation slot from the prcomp output object and plots the genes according to their absolute loading values.
#'
#' @param x Output of type prcomp. If you use mPCA, pca is a default output.
#' @param PC which PC to take the loadings from ('PCx')
#' @param n numbers of gene names to plot
#' @param return Boolean to determine whether to return a ggplot object or not (default = FALSE just returns the plot)
#'
#' @importFrom dplyr arrange mutate left_join desc
#' @importFrom tibble rownames_to_column
#' @importFrom ggplot2 ggplot geom_point annotate xlim labs
#' @importFrom cowplot theme_cowplot
#'
#' @export
#'
#'

loadingPlot <- function(x,
                        PC = "PC1",
                        n = 10,
                        return = FALSE) {
  if (class(x) != "prcomp") {
    stop("x has to be of class 'prcomp'")
  }

  if (!(class(PC) == "character" & length(PC) == 1)) {
    stop("PC must be a character vector of length 1")
  }

  gnames <- modplots::gnames

  pc_loadings <- as.data.frame(x$rotation)
  pc_loadings <- arrange(pc_loadings, desc(abs(!!as.name(PC))))
  pc_loadings <- rownames_to_column(pc_loadings, "Gene.stable.ID")
  pc_loadings <- mutate(pc_loadings, Index = 1:nrow(pc_loadings))
  pc_loadings <-
    left_join(pc_loadings, gnames, by = "Gene.stable.ID")

  ylims <- range(abs(pc_loadings[, PC]))

  lplot <- ggplot(pc_loadings,
                  aes(
                    x = Index,
                    y = abs(!!as.name(PC)),
                    label = .data$Gene.name
                  )) +
    geom_point() +
    annotate(
      geom = "text",
      x = nrow(pc_loadings) + 50,
      y = seq(ylims[2], 0.04, length.out = n),
      hjust = 1,
      label = c(pc_loadings[1:n, "Gene.name"])
    ) +
    xlim(0, nrow(pc_loadings) + 60) +
    theme_cowplot() +
    labs(y = paste0("absolute loading of ", PC))

  if (return) {
    return(lplot)
  }

  lplot

}
