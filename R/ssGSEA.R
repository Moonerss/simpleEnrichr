#' Generate single-sample gene-set enrichment score
#'
#' @description Estimate gene-set enrichment score across all samples.
#'
#' @param exp Numeric matrix containing the expression data or gene expression signatures, with samples in columns and genes in rows.
#' @param gene.list Gene sets provided either as a list object or as a GeneSetCollection object.
#' @param method Method to employ in the estimation of gene-set enrichment scores per sample.
#' By default this is set to gsva (Hänzelmann et al, 2013) and other options are ssgsea (Barbie et al, 2009), zscore (Lee et al, 2008) or plage (Tomfohr et al, 2005).
#' The latter two standardize first expression profiles into z-scores over the samples and,
#' in the case of zscore, it combines them together as their sum divided by the square-root of the size of the gene set,
#' while in the case of plage they are used to calculate the singular value decomposition (SVD) over the genes in the
#' gene set and use the coefficients of the first right-singular vector as pathway activity profile.
#' @param verbose Gives information about each calculation step. Default: \code{TRUE}.
#' @param BPPARAM An object of class [`BiocParallelParam`] specifying parameters
#'   related to the parallel execution of some of the tasks and calculations
#'   within this function.
#'
#' @importFrom GSVA gsva ssgseaParam gsvaParam plageParam plageParam
#' @importFrom BiocParallel SerialParam
#'
#' @export
#'
simple_ssgsea <- function(exp, gene.list, method = c("ssgsea", "gsva", "zscore", "plage"), verbose = TRUE, BPPARAM = SerialParam(progressbar = verbose)) {

  method <- match.arg(method)
  method_cmd <- switch(method,
                       ssgsea = ssgseaParam,
                       gsva = gsvaParam,
                       zscore = plageParam,
                       plage = plageParam
  )
  if (method == 'gsva') {
    param <- method_cmd(exprData = as.matrix(exp), geneSets = gene.list, kcdf = "Gaussian")
  } else {
    param <- method_cmd(exprData = as.matrix(exp), geneSets = gene.list)
  }
  cli::cli_alert_info('Run gene set score use method: {.val {method}}')
  score <- as.data.frame(t(GSVA::gsva(param = param, verbose = verbose, BPPARAM = BPPARAM)))

  return(score)

}
