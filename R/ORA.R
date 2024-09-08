#' Over-representative analysis
#' @description
#' Perform over-representative analysis included GO, KEGG, WikiPathways, Reactome, MsigDB, Disease Ontoloty, Cancer Gene Network, DisGeNET, CellMarker, and CMAP.
#'
#' @param genes A vector of gene id.
#' @param background.genes Background genes. If missing, the all genes listed in the database (eg TERM2GENE table) will be used as background.
#' @param gene.type gene id type input gene.
#' @param enrich.type Select an enrichment method. One of GO, KEGG, MKEGG, WikiPathways, Reactome, MsigDB, DO, CGN, DisGeNET, CellMarker, and CMAP.
#' @param organism Specify species, currently support only Human and Mouse.
#' @param GO.ont GO parameter. One of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
#' @param GO.simplify GO parameter. Whether to remove redundancy of enriched GO terms. If TRUE, will return a list with raw result and result after remove redundancy.
#' @param KEGG.use.internal.data KEGG parameter. Logical, use KEGG.db or latest online KEGG data.
#' @param MsigDB.category MsigDB parameter. MSigDB collection abbreviation, such as All, H, C1, C2, C3, C4, C5, C6, C7.
#' @param CMAP.min.Geneset.Size CMAP parameter. Minimal size of CMAP genes annotated for testing. Recommended use 3.
#' @param pvalue.cutoff pvalue cutoff on enrichment tests to report as significant.
#' @param qvalue.cutoff qvalue cutoff on enrichment tests to report as significant.
#' @param padjust.method one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' @param min.Geneset.Size Minimal size of genes annotated for testing. Not suitable for CMAP.
#' @param max.Geneset.Size Maximal size of genes annotated for testing.
#'
#' @importFrom rlang .data
#' @importFrom DOSE enrichDO enrichNCG enrichDGN
#' @importFrom ReactomePA enrichPathway
#' @importFrom msigdbr msigdbr
#' @importFrom tidyr unnest
#' @importFrom vroom vroom
#'
#' @examples
#' genes <- c("CANX", "HSPA1B", "KLRC2", "PSMC6", "RFXAP", "TAP1")
#' res <- simple_ORA(genes, enrich.type = "GO")
#'
#' @export
#'
simple_ORA <- function(genes, background.genes = NULL, gene.type = "SYMBOL",
                       enrich.type = c(
                         "GO", "KEGG", "MKEGG", "WikiPathways", "Reactome",
                         "MsigDB", "DO", "CGN", "DisGeNET", "CellMarker", "CMAP"
                       ),
                       organism = c("Human", "Mouse"),
                       GO.ont = c("BP", "CC", "MF", "ALL"), GO.simplify = T,
                       KEGG.use.internal.data = F, MsigDB.category = "H", CMAP.min.Geneset.Size = 3,
                       pvalue.cutoff = 0.05, qvalue.cutoff = 0.05, padjust.method = "BH",
                       min.Geneset.Size = 10, max.Geneset.Size = 1000) {
  ## check argument
  organism <- match.arg(organism)
  GO.ont <- match.arg(GO.ont)
  enrich.type <- match.arg(enrich.type)

  if (gene.type == "SYMBOL") {
    cli::cli_alert_info("Updating gene symbols...")
    genes <- suppressWarnings(update_symbol(genes = genes, species = organism, unmapGene_keep = T)[[2]])
    if (!is.null(background.genes)) {
      background.genes <- suppressWarnings(update_symbol(genes = stats::na.omit(background.genes), unmapGene_keep = T)[[2]])
    }
  }
  if (organism == "Human") {
    OrgDb <- "org.Hs.eg.db"
  }
  if (organism == "Mouse") {
    OrgDb <- "org.Mm.eg.db"
  }
  if (gene.type != "ENTREZID") {
    cli::cli_alert_info("Transforming {.val {gene.type}} to ENTREZID...")
    genes <- suppressWarnings(clusterProfiler::bitr(gene = genes, fromType = gene.type, toType = "ENTREZID", OrgDb = OrgDb))[, "ENTREZID"]
    if (!is.null(background.genes)) {
      background.genes <- suppressWarnings(clusterProfiler::bitr(gene = background.genes, fromType = gene.type, toType = "ENTREZID", OrgDb = OrgDb))[, "ENTREZID"]
    }
  }
  if (enrich.type == "GO") {
    cli::cli_alert_info("Performing GO-{.val {GO.ont}} enrichment...")
    ego <- suppressWarnings(suppressMessages(clusterProfiler::enrichGO(
      gene = genes,
      OrgDb = OrgDb, keyType = "ENTREZID", ont = GO.ont,
      pvalueCutoff = pvalue.cutoff, qvalueCutoff = qvalue.cutoff,
      pAdjustMethod = padjust.method, minGSSize = min.Geneset.Size,
      maxGSSize = max.Geneset.Size, universe = background.genes,
      readable = T
    )))
    if (!is.null(ego)) {
      ego <- getEF(ego)
    }
    cli::cli_alert_info("{.val {nrow(ego)}} significant terms were detected...")
    if (GO.simplify & !is.null(ego)) {
      cli::cli_alert_info("Symplifying GO results...")
      simple_ego <- clusterProfiler::simplify(ego)
      if (!is.null(simple_ego)) {
        simple_ego <- getEF(simple_ego)
      }
      cli::cli_alert_info("Return a list with raw enrich result and symplify result ...")
      res <- list(GO = ego, simplyGO = simple_ego)
    } else {
      res <- ego
    }
  }
  if (enrich.type == "KEGG") {
    cli::cli_alert_info("Performing KEGG enrichment...")
    if (organism == "Human") {
      KEGG.organism <- "hsa"
    }
    if (organism == "Mouse") {
      KEGG.organism <- "mmu"
    }
    ekegg <- suppressWarnings(suppressMessages(clusterProfiler::enrichKEGG(
      gene = genes,
      organism = KEGG.organism, pvalueCutoff = pvalue.cutoff,
      qvalueCutoff = qvalue.cutoff, pAdjustMethod = padjust.method,
      minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size,
      universe = background.genes, use_internal_data = KEGG.use.internal.data
    )))
    res <- ekegg
    if (!is.null(res)) {
      res <- getEF(res)
    }
    cli::cli_alert_info("{.val {nrow(res)}} significant terms were detected...")
  }
  if (enrich.type == "MKEGG") {
    cli::cli_alert_info("Performing Module KEGG enrichment...")
    if (organism == "Human") {
      KEGG.organism <- "hsa"
    }
    if (organism == "Mouse") {
      KEGG.organism <- "mmu"
    }
    emkegg <- suppressWarnings(suppressMessages(clusterProfiler::enrichMKEGG(
      gene = genes,
      organism = KEGG.organism, pvalueCutoff = pvalue.cutoff,
      qvalueCutoff = qvalue.cutoff, pAdjustMethod = padjust.method,
      minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size,
      universe = background.genes
    )))
    res <- emkegg
    if (!is.null(res)) {
      res <- getEF(res)
    }
    cli::cli_alert_info("{.val {nrow(res)}} significant terms were detected...")
  }
  if (enrich.type == "WikiPathways") {
    cli::cli_alert_info("Performing WikiPathways enrichment...")
    if (organism == "Human") {
      WikiPathways.organism <- "Homo sapiens"
    }
    if (organism == "Mouse") {
      WikiPathways.organism <- "Mus musculus"
    }
    eWP <- suppressWarnings(clusterProfiler::enrichWP(
      gene = genes,
      organism = WikiPathways.organism, pvalueCutoff = pvalue.cutoff,
      qvalueCutoff = qvalue.cutoff, pAdjustMethod = padjust.method,
      minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size,
      universe = background.genes
    ))
    res <- eWP
    if (!is.null(res)) {
      res <- getEF(res)
    }
    cli::cli_alert_info("{.val {nrow(res)}} significant terms were detected...")
  }
  if (enrich.type == "Reactome") {
    cli::cli_alert_info("Performing Reactome pathways enrichment...")
    if (organism == "Human") {
      Reactome.organism <- "human"
    }
    if (organism == "Mouse") {
      Reactome.organism <- "mouse"
    }
    eRP <- suppressWarnings(ReactomePA::enrichPathway(
      gene = genes,
      organism = Reactome.organism, pvalueCutoff = pvalue.cutoff,
      qvalueCutoff = qvalue.cutoff, pAdjustMethod = padjust.method,
      minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size,
      readable = T, universe = background.genes
    ))
    res <- eRP
    if (!is.null(res)) {
      res <- getEF(res)
    }
    cli::cli_alert_info("{.val {nrow(res)}} significant terms were detected...")
  }
  if (enrich.type == "DO") {
    cli::cli_alert_info("Performing Disease Ontoloty enrichment...")
    eDO <- suppressWarnings(DOSE::enrichDO(
      gene = genes,
      ont = "DO", pvalueCutoff = pvalue.cutoff, qvalueCutoff = qvalue.cutoff,
      pAdjustMethod = padjust.method, minGSSize = min.Geneset.Size,
      maxGSSize = max.Geneset.Size, readable = T, universe = background.genes
    ))
    res <- eDO
    if (!is.null(res)) {
      res <- getEF(res)
    }
    cli::cli_alert_info("{.val {nrow(res)}} significant terms were detected...")
  }
  if (enrich.type == "CGN") {
    cli::cli_alert_info("Performing Cancer Gene Network enrichment...")
    eNCG <- suppressWarnings(DOSE::enrichNCG(
      gene = genes,
      pvalueCutoff = pvalue.cutoff, qvalueCutoff = qvalue.cutoff,
      pAdjustMethod = padjust.method, minGSSize = min.Geneset.Size,
      maxGSSize = max.Geneset.Size, readable = T, universe = background.genes
    ))
    res <- eNCG
    if (!is.null(res)) {
      res <- getEF(res)
    }
    cli::cli_alert_info("{.val {nrow(res)}} significant terms were detected...")
  }
  if (enrich.type == "DisGeNET") {
    cli::cli_alert_info("Performing DisGeNET enrichment...")
    eDGN <- suppressWarnings(DOSE::enrichDGN(
      gene = genes,
      pvalueCutoff = pvalue.cutoff, qvalueCutoff = qvalue.cutoff,
      pAdjustMethod = padjust.method, minGSSize = min.Geneset.Size,
      maxGSSize = max.Geneset.Size, readable = T, universe = background.genes
    ))
    res <- eDGN
    if (!is.null(res)) {
      res <- getEF(res)
    }
    cli::cli_alert_info("{.val {nrow(res)}} significant terms were detected...")
  }
  if (enrich.type == "CellMarker") {
    cli::cli_alert_info("Performing CellMarker enrichment...")
    if (organism == "Human") {
      cell_marker_data <- suppressMessages(vroom::vroom("http://xteam.xbio.top/CellMarker/download/Human_cell_markers.txt"))
    }
    if (organism == "Mouse") {
      cell_marker_data <- suppressMessages(vroom::vroom("http://xteam.xbio.top/CellMarker/download/Mouse_cell_markers.txt"))
    }
    cells <- cell_marker_data %>%
      dplyr::select(
        .data$cellName,
        .data$geneID
      ) %>%
      dplyr::mutate(geneID = strsplit(
        .data$geneID,
        ", "
      )) %>%
      tidyr::unnest(cols = .data$geneID)
    eCM <- suppressWarnings(clusterProfiler::enricher(
      gene = genes,
      pvalueCutoff = pvalue.cutoff, qvalueCutoff = qvalue.cutoff,
      pAdjustMethod = padjust.method, minGSSize = min.Geneset.Size,
      maxGSSize = max.Geneset.Size, TERM2GENE = cells,
      universe = background.genes
    ))
    res <- eCM
    if (!is.null(res)) {
      res <- getEF(res)
    }
    cli::cli_alert_info("{.val {nrow(res)}} significant terms were detected...")
  }
  if (enrich.type == "MsigDB") {
    cli::cli_alert_info("Performing MsigDB-{.val {MsigDB.category}} enrichment...")
    if (organism == "Human") {
      MsigDB.organism <- "Homo sapiens"
    }
    if (organism == "Mouse") {
      MsigDB.organism <- "Mus musculus"
    }
    mall <- msigdbr::msigdbr(species = MsigDB.organism)
    if (MsigDB.category == "All" | MsigDB.category == "ALL") {
      mg <- mall %>% dplyr::select(.data$gs_name, .data$entrez_gene)
    } else {
      mg <- msigdbr::msigdbr(
        species = MsigDB.organism,
        category = MsigDB.category
      ) %>% dplyr::select(
        .data$gs_name,
        .data$entrez_gene
      )
    }
    eMSIG <- suppressWarnings(clusterProfiler::enricher(
      gene = genes,
      pvalueCutoff = pvalue.cutoff, qvalueCutoff = qvalue.cutoff,
      pAdjustMethod = padjust.method, minGSSize = min.Geneset.Size,
      maxGSSize = max.Geneset.Size, TERM2GENE = mg, universe = background.genes
    ))
    res <- eMSIG
    if (!is.null(res)) {
      res <- getEF(res)
    }
    cli::cli_alert_info("{.val {nrow(res)}} significant terms were detected...")
  }
  if (enrich.type == "CMAP") {
    if (organism != "Human") {
      cli::cli_abort("CMAP only supports organism = Human!")
    }
    cli::cli_alert_info("Performing CMAP enrichment...")
    eCMAP <- suppressWarnings(clusterProfiler::enricher(
      gene = genes,
      pvalueCutoff = pvalue.cutoff, qvalueCutoff = qvalue.cutoff,
      pAdjustMethod = padjust.method, minGSSize = CMAP.min.Geneset.Size,
      maxGSSize = max.Geneset.Size, TERM2GENE = CMAPfromDSEATM[
        ,
        seq_len(2)
      ], universe = background.genes
    ))
    res <- eCMAP
    if (!is.null(res)) {
      res <- getEF(res)
    }
    cli::cli_alert_info("{.val {nrow(res)}} significant terms were detected...")
  }
  if (enrich.type != "GO") {
    res <- DOSE::setReadable(res, OrgDb = OrgDb, keyType = "ENTREZID")
  }
  cli::cli_alert_success("Done!")
  return(res)
}


#' Run over-representative analysis for all gene set
#'
#' @description Perform integrated over-representative enrichment analysis included GO, KEGG, WikiPathways, Reactome, MsigDB, Disease Ontoloty, Cancer Gene Network, DisGeNET, CellMarker, and CMAP.
#' @param genes The same as in `simple_ORA`.
#' @param background.genes The same as in `simple_ORA`.
#' @param gene.type The same as in `simple_ORA`.
#' @param organism The same as in `simple_ORA`.
#' @param enrich.type A vector of selected enrichment method. You can set one or more of GO, KEGG, MKEGG, WikiPathways, Reactome, MsigDB, DO, CGN, DisGeNET, CellMarker, and CMAP.
#' @param GO.ont The same as in `simple_ORA`.
#' @param KEGG.use.internal.data The same as in `simple_ORA`.
#' @param MsigDB.category The same as in `simple_ORA`.
#' @param pvalue.cutoff The same as in `simple_ORA`.
#' @param qvalue.cutoff The same as in `simple_ORA`.
#' @param padjust.method The same as in `simple_ORA`.
#' @param min.Geneset.Size The same as in `simple_ORA`.
#' @param max.Geneset.Size The same as in `simple_ORA`.
#' @param CMAP.min.Geneset.Size The same as in `simple_ORA`.
#'
#' @examples
#' \dontrun{
#'  genes <- c("CANX", "HSPA1B", "KLRC2", "PSMC6", "RFXAP", "TAP1")
#'  res <- ORA_intergated(genes, enrich.type = c("GO", "KEGG"))
#' }
#'
#' @export
#'

ORA_intergated <- function (genes, background.genes = NULL,
                            gene.type = "SYMBOL", organism = c("Human", "Mouse"),
                            enrich.type = c("GO", "KEGG", "MKEGG", "WikiPathways", "Reactome",
                                            "MsigDB", "DO", "CGN", "DisGeNET", "CellMarker", "CMAP"),
                            GO.ont = "BP", KEGG.use.internal.data = F, MsigDB.category = "H",
                            pvalue.cutoff = 0.05, qvalue.cutoff = 0.05, padjust.method = "BH",
                            min.Geneset.Size = 10, max.Geneset.Size = 1000, CMAP.min.Geneset.Size = 3)
{

  # check method
  a <- sapply(enrich.type, function(x) {is.element(x, enrich_methods())})
  if (!all(a == TRUE)) {
    cli::cli_abort('The element in {.var enrich.type} is not all supported, please check ...')
  }

  enrich_res <- purrr::map(enrich.type, function(x) {
    simple_ORA(genes = genes, background.genes = background.genes, gene.type = gene.type,
               enrich.type = x, organism = organism,
               GO.ont = GO.ont, GO.simplify = T,
               KEGG.use.internal.data = KEGG.use.internal.data,
               MsigDB.category = MsigDB.category,
               pvalue.cutoff = pvalue.cutoff,
               qvalue.cutoff = qvalue.cutoff,
               padjust.method = padjust.method,
               min.Geneset.Size = min.Geneset.Size,
               max.Geneset.Size = max.Geneset.Size,
               CMAP.min.Geneset.Size = CMAP.min.Geneset.Size
    )
  })
  names(enrich_res) <- enrich.type
  res <- list()
  for(i in enrich.type) {
    if (i == 'GO') {
      res[['GO']] <- enrich_res$GO$GO
      res[['simplyGO']] <- enrich_res$GO$simplyGO
    } else if (i == 'KEGG') {
      res[['KEGG']] <- enrich_res$KEGG
    } else if (i == 'MKEGG') {
      res[['MKEGG']] <- enrich_res$MKEGG
    } else if (i == 'WikiPathways') {
      res[['WikiPathways']] <- enrich_res$WikiPathways
    } else if (i == 'Reactome') {
      res[['Reactome']] <- enrich_res$Reactome
    } else if (i == 'MsigDB') {
      res[['MsigDB']] <- enrich_res$MsigDB
    } else if (i == 'DO') {
      res[['DO']] <- enrich_res$DO
    } else if (i == 'CGN') {
      res[['CGN']] <- enrich_res$CGN
    } else if (i == 'DisGeNET') {
      res[['DisGeNET']] <- enrich_res$DisGeNET
    } else if (i == 'CellMarker') {
      res[['CellMarker']] <- enrich_res$CellMarker
    } else if (i == 'CMAP') {
      res[['CMAP']] <- enrich_res$CMAP
    }
  }
  cli::cli_alert_info('Final statistics ...')

  nsig <- ifelse(is.null(nrow(res$GO)), 0, nrow(res$GO)) +
    ifelse(is.null(nrow(res$KEGG)), 0, nrow(res$KEGG)) +
    ifelse(is.null(nrow(res$MKEGG)), 0, nrow(res$MKEGG)) +
    ifelse(is.null(nrow(res$WikiPathways)), 0, nrow(res$WikiPathways)) +
    ifelse(is.null(nrow(res$Reactome)), 0, nrow(res$Reactome)) +
    ifelse(is.null(nrow(res$DO)), 0, nrow(res$DO)) +
    ifelse(is.null(nrow(res$CGN)), 0, nrow(res$CGN)) +
    ifelse(is.null(nrow(res$DisGeNET)), 0, nrow(res$DisGeNET)) +
    ifelse(is.null(nrow(res$CellMarker)), 0, nrow(res$CellMarker)) +
    ifelse(is.null(nrow(res$MsigDB)), 0, nrow(res$MsigDB)) +
    ifelse(is.null(nrow(res$CMAP)), 0, nrow(res$CMAP))
  cli::cli_alert_info('{.val {nsig}} significant terms were detected...')
  cli::cli_alert_success("Done!")
  return(res)
}
