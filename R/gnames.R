#' Gallus gallus gene ID and names table (gg6)
#'
#' A data frame of chicken genes and IDs to annotate plots. Created as following:
#'
#' \preformatted{
#' gnames <- read.delim("~/Neuraltube/gg6_genes.txt")
#'    i <- sapply(gnames, is.factor)
#'    gnames[i] <- lapply(gnames[i], as.character)
#'    for (i in 1:dim(gnames)[[1]]) {
#'      ifelse(gnames[i,"Gene.name"] == "",
#'             gnames[i,"Gene.name"] <- gnames[i,"Gene.stable.ID"],
#'             gnames[i,"Gene.name"] <- gnames[i,"Gene.name"])
#'    }
#'    gnames$Gene.name <- make.unique(gnames$Gene.name, sep = "v")
#'    gnames <- gnames %>%
#'        mutate(Name = str_replace(Name, "^ENSGALG0+","EG"))
#' }
#' @format A data frame with 24356 rows and 2 variables:
#' \describe{
#'   \item{Gene.stable.ID}{Ensembl ID for chicken genes}
#'   \item{Gene.name}{corresponding Names, ID if no given name}
#'   ...
#' }
#' @source \url{http://www.https://www.ensembl.org/index.html}
"gnames"
