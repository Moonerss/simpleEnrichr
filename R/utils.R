utils::globalVariables('CMAPfromDSEATM')

# enrich methods
enrich_methods <- function() {
  c("GO", "KEGG", "MKEGG", "WikiPathways", "Reactome",
    "MsigDB", "DO", "CGN", "DisGeNET", "CellMarker", "CMAP")
}

#' Read in gene set information from .gmt files
#'
#' @description
#' This function reads in and parses information from the MSigDB's .gmt files. Pathway information will be returned as a list of gene sets.
#'
#' @param file The .gmt file to be read
#'
#' @details
#' The .gmt format is a tab-delimited list of gene sets, where each line is a separate gene set. The first column must specify the name of the gene set, and the second column is used for a short description (which this function discards). For complete details on the .gmt format, refer to the Broad Institute's Data Format's page (url: http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats).
#'
#' @return A list, where each index represents a separate gene set.
#'
#' @export
#'
read_gmt <- function(file) {

  if (!grepl("\\.gmt$", file)[1]) {
    stop("Pathway information must be a .gmt file")
  }
  geneSetDB = readLines(file)
  geneSetDB = strsplit(geneSetDB, "\t")
  names(geneSetDB) = sapply(geneSetDB, "[", 1)
  geneSetDB = lapply(geneSetDB, "[", -1:-2)
  geneSetDB = lapply(geneSetDB, function(x) {
    x[which(x != "")]
  })
  return(geneSetDB)
}




#' Calculate enrichment factor from enrichResult
#'
#' @description
#' Calculate enrichment factor from enrichResult (clusterProfiler).
#'
#' @param enrichResult enrichResult from clusterProfiler.
#'
#' @return A new enrichResult with enrichment factor.
#'
#' @examples
#' \dontrun{
#'   genes <- c("CANX", "HSPA1B", "KLRC2", "PSMC6", "RFXAP", "TAP1")
#'   obj <- clusterProfiler::enrichGO(genes, "org.Hs.eg.db",
#'                                    keyType = "SYMBOL", ont = "BP"
#'   )
#'   obj2 <- getEF(enrichResult = obj)
#'   obj2@result$EnrichFactor
#' }
#'
#' @export
#'
getEF <- function(enrichResult) {
  enrichResult@result$EnrichFactor <- apply(enrichResult@result, 1, function(x) {
    GeneRatio <- eval(parse(text = x["GeneRatio"]))
    BgRatio <- eval(parse(text = x["BgRatio"]))
    EF <- round(GeneRatio/BgRatio, 2)
  })
  return(enrichResult)
}

# export from enrichplot:::gseaScores()
gseaScores <- function(geneList, geneSet, exponent = 1, fortify = FALSE) {
  geneSet <- intersect(geneSet, names(geneList))
  N <- length(geneList)
  Nh <- length(geneSet)
  Phit <- Pmiss <- numeric(N)
  hits <- names(geneList) %in% geneSet
  Phit[hits] <- abs(geneList[hits])^exponent
  NR <- sum(Phit)
  Phit <- cumsum(Phit/NR)
  Pmiss[!hits] <- 1/(N - Nh)
  Pmiss <- cumsum(Pmiss)
  runningES <- Phit - Pmiss
  max.ES <- max(runningES)
  min.ES <- min(runningES)
  if (abs(max.ES) > abs(min.ES)) {
    ES <- max.ES
  }
  else {
    ES <- min.ES
  }
  df <- data.frame(x = seq_along(runningES), runningScore = runningES,
                   position = as.integer(hits))
  if (fortify == TRUE) {
    return(df)
  }
  df$gene = names(geneList)
  res <- list(ES = ES, runningES = df)
  return(res)
}


# export from enrichplot:::gsInfo()
gsInfo <- function (object, geneSetID) {
  geneList <- object@geneList
  if (is.numeric(geneSetID))
    geneSetID <- object@result[geneSetID, "ID"]
  geneSet <- object@geneSets[[geneSetID]]
  exponent <- object@params[["exponent"]]
  df <- gseaScores(geneList, geneSet, exponent, fortify = TRUE)
  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore))/20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList
  if (length(object@gene2Symbol) == 0) {
    df$gene <- names(geneList)
  }
  else {
    df$gene <- object@gene2Symbol[names(geneList)]
  }
  df$Description <- object@result[geneSetID, "Description"]
  return(df)
}
