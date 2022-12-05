#' Mapping Genes to enriched GO terms
#'
#' @param GOtable data.frame, goana output (filtered by topGO)
#' @param GOids ENTREZ IDs used in goana
#'
#' @importFrom org.Gg.eg.db org.Gg.egENSEMBL2EG org.Gg.egGO2ALLEGS
#' @importFrom AnnotationDbi Rkeys<- mappedLkeys Lkeys<-
#' @importFrom BiocGenerics toTable
#' @importFrom dplyr select filter rename left_join full_join
#' @importFrom tibble rownames_to_column
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @return data.frame, containing GO IDs, GOTerms, the GOtable data, and gene IDs and names.
#' @export
#'


GOgenes <- function(
    GOtable,
    GOids
    ) {

  gnames <- modplots::gnames
  # ENSEMBL to ENTREZ ID bimap
  gene_to_id <- org.Gg.egENSEMBL2EG
  # GOTerm to ENTREZ ID bimap
  go_to_id <- org.Gg.egGO2ALLEGS

  # select GO terms for IDs of genes in goana output table
  Rkeys(go_to_id) <- rownames(GOtable)
  EG <- mappedLkeys(go_to_id)
  # Genes belonging to mapped GO terms
  Lkeys(gene_to_id) <- EG

  gene_to_id <- toTable(gene_to_id)

  # create table, remove duplicates
  go_to_id <- toTable(go_to_id) %>%
    filter(.data$gene_id %in% GOids) %>%
    select(!.data$Evidence) %>%
    filter(!duplicated(.data))

  # combine GOTerms and gene ID tables, add gene names
  GO_lookup <- left_join(go_to_id, gene_to_id, by = "gene_id") %>%
    rename(Gene.stable.ID = .data$ensembl_id) %>%
    left_join(gnames, by = "Gene.stable.ID") %>%
    filter(!.data$Ontology == "CC") %>%
    filter(.data$go_id %in% rownames(GOtable))

  # combine GO table with mapped Terms and IDs
  x <- GOtable %>%
    rownames_to_column("go_id") %>%
    full_join(GO_lookup, by = "go_id") %>%
    filter(!duplicated(.data))

  x
}
