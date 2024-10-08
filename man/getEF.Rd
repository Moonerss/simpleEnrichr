% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{getEF}
\alias{getEF}
\title{Calculate enrichment factor from enrichResult}
\usage{
getEF(enrichResult)
}
\arguments{
\item{enrichResult}{enrichResult from clusterProfiler.}
}
\value{
A new enrichResult with enrichment factor.
}
\description{
Calculate enrichment factor from enrichResult (clusterProfiler).
}
\examples{
\dontrun{
  genes <- c("CANX", "HSPA1B", "KLRC2", "PSMC6", "RFXAP", "TAP1")
  obj <- clusterProfiler::enrichGO(genes, "org.Hs.eg.db",
                                   keyType = "SYMBOL", ont = "BP"
  )
  obj2 <- getEF(enrichResult = obj)
  obj2@result$EnrichFactor
}

}
