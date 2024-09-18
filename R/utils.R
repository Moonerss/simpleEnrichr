utils::globalVariables('CMAPfromDSEATM')

# enrich methods
enrich_methods <- function() {
  c("GO", "KEGG", "MKEGG", "WikiPathways", "Reactome",
    "MsigDB", "DO", "CGN", "DisGeNET", "CellMarker", "CMAP")
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
