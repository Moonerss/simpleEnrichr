#' Update gene symbol by HGNC
#' @description
#' Update the gene symbol by HGNChelper
#' @param genes a vector of gene symbol.
#' @param species A character vector of length 1, either "Human" (default) or "Mouse"
#' @param unmapGene_keep whether remove unmapped gene symbol.
#'
#' @importFrom rlang .data
#' @importFrom HGNChelper checkGeneSymbols
#' @importFrom stats na.omit
#'
#' @examples
#' human <- c("FN1", "TP53", "UNKNOWNGENE","7-Sep", "9/7", "1-Mar", "Oct4",
#'            "4-Oct", "OCT4-PG4", "C19ORF71", "C19orf71")
#' update_symbol(genes = human, unmapGene_keep = TRUE)
#'
#' @export
#'
update_symbol <- function(genes, species = c('Human', 'Mouse'), unmapGene_keep = F) {
  species <- match.arg(species)
  species <- ifelse(species == 'Human', 'human', 'mouse')
  ann <- suppressWarnings(checkGeneSymbols(genes, species = species, unmapped.as.na = unmapGene_keep)) %>%
    na.omit() %>%
    dplyr::select(-.data$Approved)
  colnames(ann) <- c("ID", "Symbol")
  if (sum(grepl("///", ann$Symbol)) > 0) {
    cli::cli_alert_warning(cli::col_br_yellow('There have genes mapped to more than one suggested symbol, we only keep the first one ...'))
    ann1 <- ann[!grepl("///", ann$Symbol), ]
    ann2 <- ann[grepl("///", ann$Symbol), ]
    ann2 <- lapply(ann2$Symbol, function(x) {
      l <- unlist(strsplit(x, " /// "))
      tmp <- ann2[ann2$Symbol == x, ] %>%
        dplyr::distinct(.data$ID, .data$Symbol, .keep_all = T)
      tmp2 <- lapply(l, function(d) {
        kk <- tmp
        kk$Symbol <- d
        return(kk)
      })
      tmp2 <- tmp2[[1]]
      return(tmp2)
    })
    ann2 <- Reduce(rbind, ann2) %>%
      dplyr::distinct(.data$ID, .data$Symbol, .keep_all = T)
    ann <- rbind(ann1, ann2) %>%
      dplyr::distinct(.data$ID, .data$Symbol, .keep_all = T)
  }
  return(ann)
}

