% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{ORA_barplot}
\alias{ORA_barplot}
\title{Plot barplot for enrihment result}
\usage{
ORA_barplot(
  enrich.obj,
  x = "RichFactor",
  color.by = "p.adjust",
  show.term.num = 15,
  label_format = 30,
  colors = c("white", "#126536"),
  color.title = color.by,
  bar.width = 0.6,
  add.bar.border = FALSE,
  y.label.position = "right",
  title = NULL,
  legend.position = "right",
  ggtheme = theme_classic(),
  ...
)
}
\arguments{
\item{enrich.obj}{An object from clusterProfiler.}

\item{x}{variable for x-axis, one of 'EnrichFactor', 'GeneRatio', 'pvalue', 'p.adjust', 'Count'.}

\item{color.by}{Variable that used to color enriched terms, one of 'GeneRatio', 'pvalue', 'p.adjust', 'Count'.}

\item{show.term.num}{A number or a list of terms. If it is a number, the first n terms will be displayed. If it is a list of terms, the selected terms will be displayed.}

\item{label_format}{a numeric value sets wrap length, alternatively a custom function to format axis labels. by default wraps names longer that 30 characters.}

\item{colors}{A color vector for the bars.}

\item{color.title}{Title of color annotation legend.}

\item{bar.width}{Width of bars.}

\item{add.bar.border}{Logical. Whether to add the black border of bars.}

\item{y.label.position}{Y label position. right, on or left.}

\item{title}{Title of the plot.}

\item{legend.position}{option of legend. 'none', 'right', 'left' or two numeric variables.}

\item{ggtheme}{ggtheme of plot.}

\item{...}{Other argument of \code{ggplot2::theme}.}
}
\description{
Plot barplot for enrihment result
}
\examples{
genes <- c("CANX", "HSPA1B", "KLRC2", "PSMC6", "RFXAP", "TAP1")
res <- simple_ORA(genes, enrich.type = "GO")
ORA_barplot(res$GO, y.label.position = 'on', x = 'RichFactor')

}
