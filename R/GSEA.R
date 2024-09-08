#' Over-representative analysis
#' @description
#' Perform over-representative analysis included GO, KEGG, WikiPathways, Reactome, MsigDB, Disease Ontoloty, Cancer Gene Network, DisGeNET, CellMarker, and CMAP.
#'
#' @param genes An order ranked geneList.
#' @param gene.type gene id type input gene.
#' @param enrich.type Select an enrichment method. One of GO, KEGG, MKEGG, WikiPathways, Reactome, MsigDB, DO, CGN, DisGeNET, CellMarker, and CMAP.
#' @param organism Specify species, currently support only Human and Mouse.
#' @param GO.ont GO parameter. One of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
#' @param GO.simplify GO parameter. Whether to remove redundancy of enriched GO terms. If TRUE, will return a list with raw result and result after remove redundancy.
#' @param KEGG.use.internal.data KEGG parameter. Logical, use KEGG.db or latest online KEGG data.
#' @param MsigDB.category MsigDB parameter. MSigDB collection abbreviation, such as All, H, C1, C2, C3, C4, C5, C6, C7.
#' @param CMAP.min.Geneset.Size CMAP parameter. Minimal size of CMAP genes annotated for testing. Recommended use 3.
#' @param pvalue.cutoff pvalue cutoff on enrichment tests to report as significant.
#' @param padjust.method one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#' @param min.Geneset.Size Minimal size of genes annotated for testing. Not suitable for CMAP.
#' @param max.Geneset.Size Maximal size of genes annotated for testing.
#'
#' @importFrom dplyr left_join
#'
#' @examples
#' \dontrun{
#' data(geneList, package="DOSE")
#'
#' # Set enrich.type using an enrichment analysis method mentioned above.
#' fit <- simple_GSEA(geneList, enrich.type = "GO", gene.type = "ENTREZID")
#' }
#'
#' @export
#'

simple_GSEA <- function(genes, gene.type = "SYMBOL",
                        enrich.type = c("GO", "KEGG", "MKEGG", "WikiPathways", "Reactome",
                                        "MsigDB", "DO", "CGN", "DisGeNET", "CellMarker", "CMAP"),
                        organism = c("Human", "Mouse"),
                        GO.ont = c("BP", "CC", "MF", "ALL"), GO.simplify = T,
                        KEGG.use.internal.data = F,
                        MsigDB.category = "H", CMAP.min.Geneset.Size = 3, pvalue.cutoff = 0.05,
                        padjust.method = "BH", min.Geneset.Size = 10, max.Geneset.Size = 1000) {

  organism <- match.arg(organism)
  GO.ont <- match.arg(GO.ont)
  enrich.type <- match.arg(enrich.type)

  d <- data.frame(raw_id = names(genes), value = genes)

  if (gene.type == "SYMBOL") {
    cli::cli_alert_info("Updating gene symbols...")
    d$raw_id <- suppressWarnings(update_symbol(genes = d$raw_id, species = organism, unmapGene_keep = T)[[2]])
    d <- d %>%
      dplyr::group_by(.data$raw_id) %>%
      dplyr::summarise_all(mean)
  }
  d <- d %>% dplyr::arrange(dplyr::desc(.data$value))
  if (organism == "Human") {
    OrgDb <- "org.Hs.eg.db"
  }
  if (organism == "Mouse") {
    OrgDb <- "org.Mm.eg.db"
  }
  if (gene.type != "ENTREZID") {
    cli::cli_alert_info("Transforming {.val {gene.type}} to ENTREZID...")
    d3 <- suppressWarnings(clusterProfiler::bitr(d$raw_id,
      fromType = gene.type, toType = "ENTREZID", OrgDb = OrgDb
    ))
    d3 <- d3 %>%
      dplyr::left_join(d, by = gene.type) %>%
      dplyr::arrange(dplyr::desc(.data$value))
    genes <- d3$value
    names(genes) <- d3$ENTREZID
  } else {
    genes <- d$value
    names(genes) <- d$raw_id
  }
  if (enrich.type == "GO") {
    cli::cli_alert_info("Performing GO-{.val {GO.ont}} enrichment...")
    ego <- suppressWarnings(suppressMessages(clusterProfiler::gseGO(
      geneList = genes,
      OrgDb = OrgDb, keyType = "ENTREZID", ont = GO.ont,
      pvalueCutoff = pvalue.cutoff, pAdjustMethod = padjust.method,
      minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size,
      eps = 1e-100, seed = T, verbose = F
    )))
    ego <- DOSE::setReadable(ego, OrgDb = OrgDb, keyType = "ENTREZID")
    cli::cli_alert_info("{.val {nrow(ego)}} significant terms were detected..")
    if (GO.simplify & !is.null(ego)) {
      cli::cli_alert_info("Symplifying GO results...")
      simple_ego <- clusterProfiler::simplify(ego)
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
    ekegg <- suppressWarnings(suppressMessages(clusterProfiler::gseKEGG(
      geneList = genes,
      organism = KEGG.organism, minGSSize = min.Geneset.Size,
      maxGSSize = max.Geneset.Size, pvalueCutoff = pvalue.cutoff,
      pAdjustMethod = padjust.method, eps = 1e-100, seed = T,
      verbose = F, use_internal_data = KEGG.use.internal.data
    )))
    res <- ekegg
    res <- DOSE::setReadable(res, OrgDb = OrgDb, keyType = "ENTREZID")
    cli::cli_alert_info("{.val {nrow(res)}} significant terms were detected..")
  }
  if (enrich.type == "MKEGG") {
    cli::cli_alert_info("Performing Module KEGG enrichment...")
    if (organism == "Human") {
      KEGG.organism <- "hsa"
    }
    if (organism == "Mouse") {
      KEGG.organism <- "mmu"
    }
    emkegg <- suppressWarnings(suppressMessages(clusterProfiler::gseMKEGG(
      geneList = genes,
      organism = KEGG.organism, keyType = "kegg", minGSSize = min.Geneset.Size,
      maxGSSize = max.Geneset.Size, pvalueCutoff = pvalue.cutoff,
      pAdjustMethod = padjust.method, eps = 1e-100, seed = T,
      verbose = F
    )))
    res <- emkegg
    res <- DOSE::setReadable(res, OrgDb = OrgDb, keyType = "ENTREZID")
    cli::cli_alert_info("{.val {nrow(res)}} significant terms were detected..")
  }
  if (enrich.type == "WikiPathways") {
    cli::cli_alert_info("+++ Performing WikiPathways enrichment...")
    if (organism == "Human") {
      WikiPathways.organism <- "Homo sapiens"
    }
    if (organism == "Mouse") {
      WikiPathways.organism <- "Mus musculus"
    }
    eWP <- suppressWarnings(clusterProfiler::gseWP(
      geneList = genes,
      organism = WikiPathways.organism, minGSSize = min.Geneset.Size,
      maxGSSize = max.Geneset.Size, pvalueCutoff = pvalue.cutoff,
      pAdjustMethod = padjust.method, eps = 1e-100, seed = T,
      verbose = F
    ))
    res <- eWP
    res <- DOSE::setReadable(res, OrgDb = OrgDb, keyType = "ENTREZID")
    cli::cli_alert_info("{.val {nrow(res)}} significant terms were detected..")
  }
  if (enrich.type == "Reactome") {
    cli::cli_alert_info("Performing Reactome pathways enrichment...")
    if (organism == "Human") {
      Reactome.organism <- "human"
    }
    if (organism == "Mouse") {
      Reactome.organism <- "mouse"
    }
    eRP <- suppressWarnings(ReactomePA::gsePathway(
      geneList = genes,
      organism = Reactome.organism, minGSSize = min.Geneset.Size,
      maxGSSize = max.Geneset.Size, pvalueCutoff = pvalue.cutoff,
      pAdjustMethod = padjust.method, eps = 1e-100, seed = T,
      verbose = F
    ))
    res <- eRP
    res <- DOSE::setReadable(res, OrgDb = OrgDb, keyType = "ENTREZID")
    cli::cli_alert_info("{.val {nrow(res)}} significant terms were detected..")
  }
  if (enrich.type == "DO") {
    cli::cli_alert_info("Performing Disease Ontoloty enrichment...")
    eDO <- suppressWarnings(DOSE::gseDO(
      geneList = genes,
      minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size,
      pvalueCutoff = pvalue.cutoff, pAdjustMethod = padjust.method,
      eps = 1e-100, seed = T, verbose = F
    ))
    res <- eDO
    res <- DOSE::setReadable(res, OrgDb = OrgDb, keyType = "ENTREZID")
    cli::cli_alert_info("{.val {nrow(res)}} significant terms were detected..")
  }
  if (enrich.type == "CGN") {
    cli::cli_alert_info("Performing Cancer Gene Network enrichment...")
    eNCG <- suppressWarnings(DOSE::gseNCG(
      geneList = genes,
      minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size,
      pvalueCutoff = pvalue.cutoff, pAdjustMethod = padjust.method,
      eps = 1e-100, seed = T, verbose = F
    ))
    res <- eNCG
    res <- DOSE::setReadable(res, OrgDb = OrgDb, keyType = "ENTREZID")
    cli::cli_alert_info("{.val {nrow(res)}} significant terms were detected..")
  }
  if (enrich.type == "DisGeNET") {
    cli::cli_alert_info("Performing DisGeNET enrichment...")
    eDGN <- suppressWarnings(DOSE::gseDGN(
      geneList = genes,
      minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size,
      pvalueCutoff = pvalue.cutoff, pAdjustMethod = padjust.method,
      eps = 1e-100, seed = T, verbose = F
    ))
    res <- eDGN
    res <- DOSE::setReadable(res, OrgDb = OrgDb, keyType = "ENTREZID")
    cli::cli_alert_info("{.val {nrow(res)}} significant terms were detected..")
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
    eCM <- suppressWarnings(clusterProfiler::GSEA(
      geneList = genes,
      minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size,
      pvalueCutoff = pvalue.cutoff, pAdjustMethod = padjust.method,
      eps = 1e-100, seed = T, verbose = F, TERM2GENE = cells
    ))
    res <- eCM
    res <- DOSE::setReadable(res, OrgDb = OrgDb, keyType = "ENTREZID")
    cli::cli_alert_info("{.val {nrow(res)}} significant terms were detected..")
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
    eMSIG <- suppressWarnings(clusterProfiler::GSEA(
      geneList = genes,
      minGSSize = min.Geneset.Size, maxGSSize = max.Geneset.Size,
      pvalueCutoff = pvalue.cutoff, pAdjustMethod = padjust.method,
      eps = 1e-100, seed = T, verbose = F, TERM2GENE = mg
    ))
    res <- eMSIG
    res <- DOSE::setReadable(res, OrgDb = OrgDb, keyType = "ENTREZID")
    cli::cli_alert_info("{.val {nrow(res)}} significant terms were detected..")
  }
  if (enrich.type == "CMAP") {
    if (organism != "Human") {
      cli::cli_abort("CMAP only supports organism = Human!")
    }
    cli::cli_alert_info("Performing CMAP enrichment...")
    eCMAP <- suppressWarnings(clusterProfiler::GSEA(
      geneList = genes,
      minGSSize = CMAP.min.Geneset.Size, maxGSSize = max.Geneset.Size,
      pvalueCutoff = pvalue.cutoff, pAdjustMethod = padjust.method,
      eps = 1e-100, seed = T, verbose = F, TERM2GENE = CMAPfromDSEATM[
        ,
        seq_len(2)
      ]
    ))
    res <- eCMAP
    res <- DOSE::setReadable(res, OrgDb = OrgDb, keyType = "ENTREZID")
    cli::cli_alert_info("{.val {nrow(res)}} significant terms were detected..")
  }
  cli::cli_alert_success("Done!")
  return(res)
}

#' Run GSEA analysis for all gene set
#'
#' @description Perform integrated gene set enrichment analysis included GO, KEGG, WikiPathways, Reactome, MsigDB, Disease Ontoloty, Cancer Gene Network, DisGeNET, CellMarker, and CMAP.
#' @param genes The same as in `simple_GSEA`.
#' @param gene.type The same as in `simple_GSEA`.
#' @param organism The same as in `simple_GSEA`.
#' @param enrich.type A vector of selected enrichment method. You can set one or more of GO, KEGG, MKEGG, WikiPathways, Reactome, MsigDB, DO, CGN, DisGeNET, CellMarker, and CMAP.
#' @param GO.ont The same as in `simple_GSEA`.
#' @param KEGG.use.internal.data The same as in `simple_GSEA`.
#' @param MsigDB.category The same as in `simple_GSEA`.
#' @param pvalue.cutoff The same as in `simple_GSEA`.
#' @param padjust.method The same as in `simple_GSEA`.
#' @param min.Geneset.Size The same as in `simple_GSEA`.
#' @param max.Geneset.Size The same as in `simple_GSEA`.
#' @param CMAP.min.Geneset.Size The same as in `simple_GSEA`.
#'
#' @importFrom purrr map
#'
#' @examples
#' \dontrun{
#'   data(geneList, package="DOSE")
#'
#'   # Set enrich.type using an enrichment analysis method mentioned above.
#'   fit <- GSEA_intergated(geneList, gene.type = 'ENTREZID', enrich.type = c("GO", "KEGG"))
#' }
#'
#' @export
#'

GSEA_intergated <- function (genes, gene.type = "SYMBOL", organism = c("Human", "Mouse"),
                             enrich.type = c("GO", "KEGG", "MKEGG", "WikiPathways", "Reactome",
                                             "MsigDB", "DO", "CGN", "DisGeNET", "CellMarker", "CMAP"),
                             GO.ont = "BP", KEGG.use.internal.data = F, MsigDB.category = "H",
                             pvalue.cutoff = 0.05, padjust.method = "BH",
                             min.Geneset.Size = 10, max.Geneset.Size = 1000, CMAP.min.Geneset.Size = 3)
{

  # check method
  a <- sapply(enrich.type, function(x) {is.element(x, enrich_methods())})
  if (!all(a == TRUE)) {
    cli::cli_abort('The element in {.var enrich.type} is not all supported, please check ...')
  }

  enrich_res <- purrr::map(enrich.type, function(x) {
    simple_GSEA(genes = genes, gene.type = gene.type,
                enrich.type = x, organism = organism,
                GO.ont = GO.ont, GO.simplify = T,
                KEGG.use.internal.data = KEGG.use.internal.data,
                MsigDB.category = MsigDB.category,
                pvalue.cutoff = pvalue.cutoff,
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
