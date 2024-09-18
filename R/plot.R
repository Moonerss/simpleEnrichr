#' Plot Dotplot for enrihment result
#' @description Plot Dotplot for enrihment result
#'
#' @param enrich.obj An object from clusterProfiler.
#' @param x variable for x-axis, one of 'GeneRatio', 'pvalue', 'p.adjust', 'Count'.
#' @param color.by Variable that used to color enriched terms, one of 'GeneRatio', 'pvalue', 'p.adjust', 'Count'.
#' @param show.term.num A number or a list of terms. If it is a number, the first n terms will be displayed. If it is a list of terms, the selected terms will be displayed.
#' @param label_format a numeric value sets wrap length, alternatively a custom function to format axis labels. by default wraps names longer that 30 characters.
#' @param colors A color vector for the bars.
#' @param color.title Title of color annotation legend.
#' @param size.by Variable that used to size enriched terms, one of 'GeneRatio', 'pvalue', 'p.adjust', 'Count'
#' @param size.range Two numeric variables, the first is minimal value and the first is maximal value.
#' @param size.title Title of size annotation legend.
#' @param y.label.position Y label position. right or left.
#' @param title Title of the plot.
#' @param legend.position ostion of legend. 'none', 'right', 'left' or two numeric variables.
#' @param ggtheme ggtheme of plot.
#' @param ... Other argument of `ggplot2::theme`.
#'
#' @import ggplot2
#' @importFrom dplyr desc
#' @importFrom Hmisc capitalize
#'
#' @examples
#' genes <- c("CANX", "HSPA1B", "KLRC2", "PSMC6", "RFXAP", "TAP1")
#' res <- simple_ORA(genes, enrich.type = "GO")
#' ORA_dotplot(res$simplyGO)
#'
#' @export
#'
ORA_dotplot <- function(enrich.obj,
                        x = "GeneRatio", color.by = "p.adjust",
                        show.term.num = 15, label_format = 30,
                        colors = c("blue", "white", "red"),
                        color.title = color.by, size.by = "Count",
                        size.range = c(3, 8),
                        size.title = size.by,
                        y.label.position = "right",
                        title = NULL, legend.position = "right",
                        ggtheme = theme_bw(base_rect_size = 1),
                        ...) {

  enrich.obj@result$GeneRatio <- apply(enrich.obj@result, 1, function(x) {
    eval(parse(text = x["GeneRatio"]))
  })
  x.lab <- x
  ## x axis
  if (x == "pvalue") {
    enrich.obj@result$Sig <- -log10(enrich.obj@result$pvalue)
    x.lab <- bquote(~ -Log[10] ~ italic("P-value"))
    x <- "Sig"
  } else if (x == "p.adjust") {
    enrich.obj@result$Sig <- -log10(enrich.obj@result$p.adjust)
    x.lab <- bquote(~ -Log[10] ~ "FDR")
    x <- "Sig"
  }

  ## color
  if (color.by == "pvalue") {
    enrich.obj@result$SigL <- -log10(enrich.obj@result$pvalue)
    color.title <- bquote(~ -Log[10] ~ italic("P-value"))
    color.by <- "SigL"
  }
  if (color.by == "p.adjust") {
    enrich.obj@result$SigL <- -log10(enrich.obj@result$p.adjust)
    color.title <- bquote(~ -Log[10] ~ "FDR")
    color.by <- "SigL"
  }

  ## size
  if (size.by == "pvalue") {
    enrich.obj@result$SigL2 <- -log10(enrich.obj@result$pvalue)
    size.title <- bquote(~ -Log[10] ~ italic("P-value"))
    size.by <- "SigL2"
  } else if (size.by == "p.adjust") {
    enrich.obj@result$SigL2 <- -log10(enrich.obj@result$p.adjust)
    size.title <- bquote(~ -Log[10] ~ "FDR")
    size.by <- "SigL2"
  }

  ## select term
  show.term.num <- ifelse(nrow(enrich.obj@result) >= show.term.num,
    show.term.num, nrow(enrich.obj@result)
  )

  dd <- enrich.obj@result %>%
    dplyr::arrange(.data$pvalue)
  dd <- dd[1:show.term.num,]%>%
    dplyr::arrange(desc(get(x)))

  dd <- dplyr::mutate(dd, Description = factor(.data$Description, rev(.data$Description)))
  p <- ggplot(dd) +
    aes(.data[[x]], .data[['Description']], fill = .data[[color.by]], size = .data[[size.by]]) +
    geom_point(shape = 21, color = "black") +
    labs(fill = color.title, size = size.title, x = x.lab, y = NULL, title = title) +
    scale_fill_gradientn(colours = colors) +
    scale_size(range = size.range) +
    scale_y_discrete(labels = Hmisc::capitalize, position = y.label.position) +
    ggtheme +
    theme(axis.text.x = element_text(size = 10, colour = "black"),
          axis.title.x = element_text(size = 13, colour = "black", face = "bold"),
          axis.text.y = element_text(size = 13, colour = "black"),
          axis.title = element_blank(),
          panel.background = element_blank(),
          legend.text = element_text(size = 11, colour = "black"),
          legend.title = element_text(size = 13, colour = "black", face = "bold"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.position = legend.position,
          plot.title = element_text(hjust = 0.5, size = 14, colour = "black", face = "bold"),
          ...)
  return(p)
}



#' Enrichment dotplot for positive or negative GSEA results
#' @description
#' Plot enrichment dotplot for positive or negative GSEA results.
#'
#' @param enrich.obj A GSEA enrichment object from `clusterProfiler`.
#' @param type Specify whether you want to show positive or negative results.
#' @param show.term.num A number or a list of terms. If it is a number, the first n terms will be displayed. If it is a list of terms, the selected terms will be displayed.
#' @param Selct.P Selct.P P value (P) or adjust P value (FDR) were selected to define significant terms.
#' @param cutoff.P A cutoff value for `Select.P`.
#' @param colors A color vector for the dots.
#' @param size.range Two numeric variables, the first is minimal value and the first is maximal value.
#' @param y.label.position Y label position. right or left.
#' @param title Title of the plot.
#' @param legend.position Position of legend. 'none', 'right', 'left' or two numeric variables.
#' @param ggtheme ggtheme of plot.
#' @param ... Other argument of `ggplot2::theme`.
#'
#' @import ggplot2
#' @importFrom dplyr desc
#' @importFrom Hmisc capitalize
#'
#' @examples
#' \dontrun{
#' data(geneList, package="DOSE")
#' GSEA_res <- simple_GSEA(geneList, gene.type = "ENTREZID", enrich.type = 'GO')
#' GSEA_dotplot(res$simplyGO)
#' }
#'
#'
#'
#' @export
#'

GSEA_dotplot <- function(enrich.obj,
                         type = c("pos", "neg"), show.term.num = 15,
                         Selct.P = c("FDR", "P"), cutoff.P = 0.05,
                         colors = c("blue", "white", "red"),
                         size.range = c(3, 8),
                         y.label.position = "right",
                         title = NULL, legend.position = "right",
                         ggtheme = theme_bw(base_rect_size = 1),
                         ...) {
  type <- match.arg(type)
  Selct.P <- match.arg(Selct.P)

  select_p <- ifelse(Selct.P == 'FDR', 'p.adjust', 'pvalue')

  if (type == 'pos') {
    r <- enrich.obj@result %>%
      dplyr::filter(.data[[select_p]] < cutoff.P, .data$NES > 0) %>%
      dplyr::mutate(sig = -log10(.data[[select_p]])) %>%
      dplyr::arrange(dplyr::desc(.data$sig))
  } else {
    r <- enrich.obj@result %>%
      dplyr::filter(.data[[select_p]] < cutoff.P, .data$NES < 0) %>%
      dplyr::mutate(sig = -log10(.data[[select_p]])) %>%
      dplyr::arrange(dplyr::desc(.data$sig))
  }

  p_str <- ifelse(Selct.P == 'FDR', 'FDR', 'P-value')
  x.lab <- bquote(~-Log[10] ~ italic(.(p_str)))

  ## select term
  show.term.num <- ifelse(nrow(r) >= show.term.num, show.term.num, nrow(r))

  if (show.term.num == 0) {
    cli::cli_abort('No significant {.val {type}} term was enriched in enrich.obj!')
  }

  dd <- r[1:show.term.num, ] %>%
    dplyr::mutate(Description = factor(.data$Description, rev(.data$Description)))

  color.title <- ifelse(type == 'pos', "NES", "abs(NES)")

  p <- ggplot(dd) +
    aes(.data$sig, .data$Description, fill = abs(.data$NES), size = abs(.data$NES)) +
    geom_point(shape = 21, color = "black") +
    labs(size = color.title, fill = color.title, x = x.lab, y = NULL, title = title) +
    scale_fill_gradientn(colours = colors) +
    scale_size(range = size.range) +
    scale_y_discrete(labels = Hmisc::capitalize, position = y.label.position) +
    guides(fill = guide_legend()) +
    ggtheme +
    theme(axis.text.x = element_text(size = 10,colour = "black"),
          axis.title.x = element_text(size = 13, colour = "black", face = "bold"),
          axis.text.y = element_text(size = 13, colour = "black"),
          axis.title.y = element_blank(),
          panel.background = element_blank(),
          legend.text = element_text(size = 11, colour = "black"),
          legend.title = element_text(size = 13, colour = "black", face = "bold"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.position = legend.position,
          plot.title = element_text(size = 14, colour = "black", face = "bold", hjust = 0.5),
          ...)
  return(p)
}


#' Visualize analyzing result of GSEA
#' @description
#' Visualize analyzing result of GSEA.
#'
#' @param GSEA.result GSEA results from `clusterProfiler::GSEA()` function.
#' @param Pathway.ID Pathway ID in the `ID` column of `GSEA.result` .
#' @param show.heatbar Whether show heat bar. Default TRUE.
#' @param show.rank Whether show Rank map. Default TRUE.
#' @param line.color Line color for running score.
#' @param rank.colors Color scheme of rank lines. A vector.
#' @param heatbar.colors Color scheme of heatbar. A vector.
#' @param add.x.ann Whether to add the title, text, and ticks of X axis.
#' @param x.lab x axis title.
#' @param line.y.lab Y label of running score plot.
#' @param rank.y.lab Y label of rank plot.
#' @param statistic.position Position of statistics in the running score plot.
#' @param statistic.face Font face of statistics.
#' @param statistic.size Font size of statistics.
#' @param ggtheme A theme object from ggplot2.
#'
#' @import enrichplot
#' @importFrom cowplot plot_grid
#'
#' @examples
#' data(geneList, package="DOSE")
#' GSEA_res <- simple_GSEA(geneList, gene.type = "ENTREZID", enrich.type = 'KEGG')
#' GSEA_rankplot(GSEA_res, Pathway.ID = 'hsa04110')
#'
#' @export
#'
GSEA_rankplot <- function(GSEA.result, Pathway.ID,
                         show.heatbar = T, show.rank = T,
                         line.color = "#41A98E",
                         rank.colors = c("#440154FF", "#482878FF", "#3E4A89FF", "#31688EFF", "#26828EFF",
                                         "#1F9E89FF", "#35B779FF", "#6DCD59FF", "#B4DE2CFF", "#FDE725FF"),
                         heatbar.colors = c("#08519C", "#3182BD", "#6BAED6", "#BDD7E7", "#EFF3FF",
                                            "#FEE5D9", "#FCAE91", "#FB6A4A", "#DE2D26", "#A50F15"),
                         add.x.ann = T,
                         x.lab = "Gene ranks",
                         line.y.lab = "Enrichment score",
                         rank.y.lab = "logFC",
                         statistic.position = c(0.5, 0.2),
                         statistic.face = "italic",
                         statistic.size = 3.5,
                         ggtheme = theme_bw(base_rect_size = 1.5)) {
  gsdata <- gsInfo(GSEA.result, Pathway.ID)
  if (!is.element(Pathway.ID, GSEA.result@result$ID)) {
    cli::cli_abort('The pathway {.val {Pathway.ID}} is not in {.var GSEA.result}, please check ...')
  }

  pathway_idx <- which(GSEA.result@result$ID == Pathway.ID)

  NES_val <- GSEA.result@result$NES[pathway_idx]
  p.adjust_val <- GSEA.result@result$p.adjust[pathway_idx]
  title <- GSEA.result@result$Description[pathway_idx]

  label <- paste0("NES = ", sprintf("%.3f", NES_val),
                  "\nFDR = ",
                  ifelse(p.adjust_val < 0.001,
                         format(p.adjust_val, scientific = T, digit = 3),
                         sprintf("%.3f", p.adjust_val)))

  ## plot
  plotlist <- list()

  plotlist[["line"]] <- ggplot(gsdata, aes(.data$x)) +
    geom_line(aes(y = .data$runningScore), linewidth = 1, color = line.color) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0.01)) +
    labs(x = NULL, y = line.y.lab,
         title = Hmisc::capitalize(title)) +
    annotate("text",
             x = nrow(gsdata) * statistic.position[1],
             y = min(gsdata$runningScore) + (max(gsdata$runningScore) - min(gsdata$runningScore)) * statistic.position[2],
             label = label, hjust = 0,
             fontface = statistic.face, size = statistic.size) +
    ggtheme +
    theme(plot.title = element_text(hjust = 0.5, size = 13, colour = "black", face = "bold"),
          panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 10, colour = "black"),
          axis.title.y = element_text(size = 13, colour = "black", face = "bold"),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          plot.margin = margin(t = 0.2, r = 0.2, b = -0.07, l = 0.2, unit = "cm"))
  if (show.heatbar) {
    plotlist[["heatbar"]] <- ggplot(gsdata, aes(.data$x)) +
      geom_linerange(aes(ymin = .data$ymin, ymax = .data$ymax), color = "grey30") +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0,0)) +
      labs(x = NULL, y = NULL) +
      ggtheme +
      theme(legend.position = "none",
            panel.grid = element_blank(),
            plot.margin = margin(t = -0.1, b = 0, unit = "cm"),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.line.x = element_blank())
    v <- seq(1, sum(gsdata$position), length.out = length(heatbar.colors))
    inv <- findInterval(rev(cumsum(gsdata$position)), v)
    if (min(inv) == 0) {
      inv <- inv + 1
    }
    col <- heatbar.colors
    ymin <- min(plotlist[["heatbar"]]$data$ymin)
    yy <- max(plotlist[["heatbar"]]$data$ymax - plotlist[["heatbar"]]$data$ymin) *
      0.3
    xmin <- which(!duplicated(inv))
    xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
    d <- data.frame(ymin = ymin, ymax = yy, xmin = xmin,
                    xmax = xmax, col = col[unique(inv)])
    plotlist[["heatbar"]] <- plotlist[["heatbar"]] +
      geom_rect(aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin, ymax = 0, fill = I(col)), data = d, alpha = 0.9, inherit.aes = FALSE)
  }

  if (show.rank) {
    plotlist[["rank"]] <- ggplot(gsdata, aes(.data$x)) +
      labs(y = rank.y.lab) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      geom_segment(aes(x = .data$x, xend = .data$x, y = .data$geneList, yend = 0, color = .data$geneList)) +
      scale_color_gradientn(colours = rank.colors) +
      ggtheme +
      theme(plot.title = element_text(hjust = 0.5, size = 13, colour = "black", face = "bold"),
            panel.grid = element_blank(),
            legend.position = "none",
            plot.margin = margin(t = -0.17, r = 0.2, b = 0.2, l = 0.2, unit = "cm"))
  }

  n <- length(plotlist)
  if (add.x.ann) {
    plotlist[[n]] <- plotlist[[n]] +
      labs(x = x.lab) +
      theme(axis.line.x = element_line(),
            axis.ticks.x = element_line(),
            axis.text.x = element_text(size = 10, colour = "black"),
            axis.title.x = element_text(size = 13, colour = "black", face = "bold"))
  }

  if (n == 3) {
    rel.heights <- c(1.5, 0.2, 1)
    p <- cowplot::plot_grid(plotlist = plotlist, ncol = 1, align = "v", rel_heights = rel.heights)
  }
  if (n == 2 & show.heatbar) {
    rel.heights <- c(1.5, 0.5)
    p <- cowplot::plot_grid(plotlist = plotlist, ncol = 1, align = "v", rel_heights = rel.heights)
  }
  if (n == 2 & show.rank) {
    rel.heights <- c(1.5, 1)
    p <- cowplot::plot_grid(plotlist = plotlist, ncol = 1, align = "v", rel_heights = rel.heights)
  }
  p
}

#' Enrichment barplot for two ORA enrichment objects
#' @description
#' Plot enrichment barplot for two ORA enrichment objects.
#'
#' @param enrich.obj1 An `enrichResult` object from clusterProfiler.
#' @param enrich.obj2 An `enrichResult` object from clusterProfiler.
#' @param Selct.P P value (P) or adjust P value (FDR) were selected to define significant terms.
#' @param cutoff.P A cutoff value for `Select.P`.
#' @param obj.types Two characters for defining the types of two objects.
#' @param obj.type.colors Two colors for the types of two objects.
#' @param obj1.top.pathway.num The number of top pathways in object 1. Based on the significant test.
#' @param obj2.top.pathway.num The number of top pathways in object 2. Based on the significant test.
#' @param bar.width Width of bar in the plot.
#' @param add.bar.border Logical. Whether to add the black border of bars.
#' @param x.limit.fold Specify the fold of x limitation.
#' @param label.size Fontsize of label.
#' @param legend.position none, left, right, top, bottom; Or Two numeric variables indicated x and y positions, respectively.
#'
#' @examples
#' genes <- c("CANX", "HSPA1B", "KLRC2", "PSMC6", "RFXAP", "TAP1")
#' res1 <- simple_ORA(genes, enrich.type = "KEGG")
#' res2 <- simple_ORA(genes, enrich.type = "MsigDB")
#' ORA_two_barplot(res1, res2, Selct.P = 'P', obj.types = c("KEGG", "MsigDB"))
#'
#' @export
#'
ORA_two_barplot <- function(enrich.obj1, enrich.obj2,
                            Selct.P = c("FDR", "P"), cutoff.P = 0.05,
                            obj.types = c("Up", "Down"),
                            obj.type.colors = c("#ED6355", "#3E94B5"),
                            obj1.top.pathway.num = 10, obj2.top.pathway.num = 10,
                            bar.width = 0.6, add.bar.border = T, x.limit.fold = 1.05,
                            label.size = 3.5, legend.position = "bottom") {
  Selct.P <- match.arg(Selct.P)
  select_p <- ifelse(Selct.P == 'FDR', 'p.adjust', 'pvalue')

  r1 <- enrich.obj1@result %>%
    dplyr::filter(.data[[select_p]] < cutoff.P) %>%
    dplyr::arrange(.data[[select_p]]) %>%
    dplyr::mutate(Type = obj.types[1])
  obj1.top.pathway.num <- ifelse(nrow(r1) >= obj1.top.pathway.num, obj1.top.pathway.num, nrow(r1))
  if (obj1.top.pathway.num == 0) {
    cli::cli_abort("No significant term was enriched in {.val enrich.obj1}")
  }
  r1 <- r1[1:obj1.top.pathway.num, ]
  r2 <- enrich.obj2@result %>%
    dplyr::filter(.data[[select_p]] < cutoff.P) %>%
    dplyr::arrange(.data[[select_p]]) %>%
    dplyr::mutate(Type = obj.types[2])

  obj2.top.pathway.num <- ifelse(nrow(r2) >= obj2.top.pathway.num, obj2.top.pathway.num, nrow(r2))
  if (obj2.top.pathway.num == 0) {
    cli::cli_abort("No significant term was enriched in {.val enrich.obj2}")
  }
  r2 <- r2[1:obj2.top.pathway.num, ]

  rr <- rbind(r1[, intersect(colnames(r1), colnames(r2))],
              r2[, intersect(colnames(r1), colnames(r2))]) %>%
    dplyr::mutate(Description = Hmisc::capitalize(.data$Description))
  rr <- rr %>%
    dplyr::mutate(log10P = log10(.data[[select_p]])) %>%
    dplyr::mutate(log10P = ifelse(.data$Type == obj.types[1], -.data$log10P, .data$log10P))

  p_str <- ifelse(Selct.P == 'FDR', 'FDR', 'P-value')
  x.lab <- bquote(~-Log[10] ~ italic(.(p_str)))

  rr <- rr %>%
    dplyr::mutate(Type = factor(.data$Type, obj.types)) %>%
    dplyr::arrange(desc(.data$log10P)) %>%
    dplyr::distinct(.data$Description, .keep_all = T) %>%
    dplyr::mutate(Description = factor(.data$Description, rev(.data$Description)))

  p <- ggplot(rr, aes(.data$log10P, .data$Description, fill = .data$Type)) +
    geom_bar(stat = "identity", width = bar.width, color = ifelse(add.bar.border, "black", NA)) +
    scale_x_continuous(limits = c(-max(abs(rr$log10P)) * x.limit.fold, max(abs(rr$log10P)) * x.limit.fold), labels = abs, name = x.lab) +
    theme_classic(base_line_size = 0.9) +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          axis.text.x = element_text(size = 10, colour = "black"),
          axis.title.x = element_text(size = 13, colour = "black", face = "bold"),
          legend.position = legend.position,
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 13, colour = "black"))

  p <- p + geom_text(data = dplyr::filter(rr, .data$Type == obj.types[2]), aes(x = 0, y = .data$Description, label = paste0(" ", .data$Description)), size = label.size, hjust = 0)
  p <- p + geom_text(data = dplyr::filter(rr, .data$Type == obj.types[1]), aes(x = -0.1, y = .data$Description, label = .data$Description), size = label.size, hjust = 1)
  p <- p + scale_fill_manual(values = obj.type.colors)

  return(p)
}


#' Enrichment barplot for two ORA enrichment objects
#' @description
#' Plot enrichment barplot for two ORA enrichment objects.
#'
#' @param enrich.obj An `gseaResult` object from clusterProfiler.
#' @param Selct.P P value (P) or adjust P value (FDR) were selected to define significant terms.
#' @param cutoff.P A cutoff value for `Select.P`.
#' @param type.colors Two colors for the positive and negative enrichment.
#' @param pos.top.pathway.num The number of top pathways in positive terms. Based on the significant test.
#' @param neg.top.pathway.num The number of top pathways in negative terms. Based on the significant test.
#' @param bar.width Width of bar in the plot.
#' @param add.bar.border Logical. Whether to add the black border of bars.
#' @param x.limit.fold Specify the fold of x limitation.
#' @param label.size Fontsize of label.
#' @param legend.position none, left, right, top, bottom; Or Two numeric variables indicated x and y positions, respectively.
#'
#' @examples
#' data(geneList, package="DOSE")
#'
#' # Set enrich.type using an enrichment analysis method mentioned above.
#' fit <- simple_GSEA(geneList, enrich.type = "KEGG", gene.type = "ENTREZID")
#' GSEA_two_barplot(fit)
#'
#' @export
#'
GSEA_two_barplot <- function(enrich.obj, Selct.P = c("FDR", "P"), cutoff.P = 0.05,
                             type.colors = c("#ED6355", "#3E94B5"),
                             pos.top.pathway.num = 10,
                             neg.top.pathway.num = 10,
                             bar.width = 0.6, add.bar.border = T,
                             x.limit.fold = 1.05, label.size = 3.5,
                             legend.position = "bottom") {
  Selct.P <- match.arg(Selct.P)
  select_p <- ifelse(Selct.P == 'FDR', 'p.adjust', 'pvalue')

  types <- c("Positive", "Negative")
  r1 <- enrich.obj@result %>%
    dplyr::filter(.data[[select_p]] < cutoff.P, .data$NES > 0) %>%
    dplyr::arrange(dplyr::desc(.data$NES)) %>%
    dplyr::mutate(Type = types[1])

  pos.top.pathway.num <- ifelse(nrow(r1) >= pos.top.pathway.num, pos.top.pathway.num, nrow(r1))
  if (pos.top.pathway.num == 0) {
    cli::cli_alert_warning("No significant positive term was enriched in {.val enrich.obj}")
    r1 <- NULL
  } else {
    r1 <- r1[1:neg.top.pathway.num,]
  }

  r2 <- enrich.obj@result %>%
    dplyr::filter(.data[[select_p]] < cutoff.P, .data$NES < 0) %>%
    dplyr::arrange(.data$NES) %>%
    dplyr::mutate(Type = types[2])
  neg.top.pathway.num <- ifelse(nrow(r2) >= neg.top.pathway.num, neg.top.pathway.num, nrow(r2))
  if (neg.top.pathway.num == 0) {
    cli::cli_alert_warning("No significant negative term was enriched in {.val enrich.obj}")
    r2 <- NULL
  } else {
    r2 <- r2[1:neg.top.pathway.num,]
  }

  if (is.null(r1)) {
    if (is.null(r2)) {
      cli::cli_abort('No significant term to plot')
    } else {
      rr <- r2
      type.colors <- type.colors[2]
    }
  } else {
    if (is.null(r2)) {
      rr <- r1
      type.colors <- type.colors[1]
    } else {
      rr <- rbind(r1[, intersect(colnames(r1), colnames(r2))],
                  r2[, intersect(colnames(r1), colnames(r2))]) %>%
        dplyr::mutate(Description = Hmisc::capitalize(.data$Description))
    }
  }
  rr <- rr %>%
    dplyr::mutate(Type = factor(.data$Type, types)) %>%
    dplyr::arrange(desc(.data$NES)) %>%
    dplyr::distinct(.data$Description, .keep_all = T) %>%
    dplyr::mutate(Description = factor(.data$Description, rev(.data$Description)))

  p <- ggplot(rr, aes(.data$NES, .data$Description, fill = .data$Type)) +
    geom_bar(stat = "identity", width = bar.width, color = ifelse(add.bar.border, "black", NA)) +
    scale_x_continuous(limits = c(-max(abs(rr$NES)) * x.limit.fold, max(abs(rr$NES)) * x.limit.fold), name = "NES") +
    theme_classic(base_line_size = 0.9) +
    ggplot2::theme(axis.text.y = element_blank(),
                   axis.title.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   axis.line.y = element_blank(),
                   axis.text.x = element_text(size = 10, colour = "black"),
                   axis.title.x = element_text(size = 13, colour = "black", face = "bold"),
                   legend.position = legend.position,
                   legend.background = element_blank(),
                   legend.key = element_blank(),
                   legend.title = element_blank(),
                   legend.text = element_text(size = 13, colour = "black"))

  if (nrow(dplyr::filter(rr, .data$Type == 'Positive')) > 0) {
    p <- p +
      geom_text(data = dplyr::filter(rr, .data$Type == types[1]), aes(x = -0.1, y = .data$Description, label = .data$Description), size = label.size, hjust = 1)
  }

  if (nrow(dplyr::filter(rr, .data$Type == 'Negative')) > 0) {
    p <- p +
      geom_text(data = dplyr::filter(rr, .data$Type == types[2]), aes(x = 0, y = .data$Description, label = paste0(" ", .data$Description)), size = label.size, hjust = 0)
  }

  p <- p + scale_fill_manual(values = type.colors)

  return(p)
}


#' Plot barplot for enrihment result
#' @description Plot barplot for enrihment result
#'
#' @param enrich.obj An object from clusterProfiler.
#' @param x variable for x-axis, one of 'EnrichFactor', 'GeneRatio', 'pvalue', 'p.adjust', 'Count'.
#' @param color.by Variable that used to color enriched terms, one of 'GeneRatio', 'pvalue', 'p.adjust', 'Count'.
#' @param show.term.num A number or a list of terms. If it is a number, the first n terms will be displayed. If it is a list of terms, the selected terms will be displayed.
#' @param label_format a numeric value sets wrap length, alternatively a custom function to format axis labels. by default wraps names longer that 30 characters.
#' @param colors A color vector for the bars.
#' @param color.title Title of color annotation legend.
#' @param bar.width Width of bars.
#' @param add.bar.border Logical. Whether to add the black border of bars.
#' @param y.label.position Y label position. right, on or left.
#' @param title Title of the plot.
#' @param legend.position option of legend. 'none', 'right', 'left' or two numeric variables.
#' @param ggtheme ggtheme of plot.
#' @param ... Other argument of `ggplot2::theme`.
#'
#' @import ggplot2
#' @importFrom dplyr desc
#' @importFrom Hmisc capitalize
#'
#' @examples
#' genes <- c("CANX", "HSPA1B", "KLRC2", "PSMC6", "RFXAP", "TAP1")
#' res <- simple_ORA(genes, enrich.type = "GO")
#' ORA_barplot(res$GO, y.label.position = 'on', x = 'RichFactor')
#'
#' @export
#'
ORA_barplot <- function(enrich.obj,
                        x = "RichFactor", color.by = "p.adjust",
                        show.term.num = 15, label_format = 30,
                        colors = c('white', '#126536'),
                        color.title = color.by,
                        bar.width = 0.6, add.bar.border = FALSE,
                        y.label.position = "right",
                        title = NULL, legend.position = "right",
                        ggtheme = theme_classic(),
                        ...) {

  # calculate enrich factor
  if (!is.element('EnrichFactor',colnames(enrich.obj@result))) {
    enrich.obj <- getEF(enrich.obj)
  }

  enrich.obj@result$GeneRatio <- apply(enrich.obj@result, 1, function(x) {
    eval(parse(text = x["GeneRatio"]))
  })
  x.lab <- x
  ## x axis
  if (x == "pvalue") {
    enrich.obj@result$Sig <- -log10(enrich.obj@result$pvalue)
    x.lab <- bquote(~ -Log[10] ~ italic("P-value"))
    x <- "Sig"
  } else if (x == "p.adjust") {
    enrich.obj@result$Sig <- -log10(enrich.obj@result$p.adjust)
    x.lab <- bquote(~ -Log[10] ~ "FDR")
    x <- "Sig"
  }

  ## color
  if (color.by == "pvalue") {
    enrich.obj@result$SigL <- -log10(enrich.obj@result$pvalue)
    color.title <- bquote(~ -Log[10] ~ italic("P-value"))
    color.by <- "SigL"
  }
  if (color.by == "p.adjust") {
    enrich.obj@result$SigL <- -log10(enrich.obj@result$p.adjust)
    color.title <- bquote(~ -Log[10] ~ "FDR")
    color.by <- "SigL"
  }


  ## select term
  show.term.num <- ifelse(nrow(enrich.obj@result) >= show.term.num,
                          show.term.num, nrow(enrich.obj@result)
  )

  dd <- enrich.obj@result %>%
    dplyr::arrange(.data$pvalue)
  dd <- dd[1:show.term.num,]%>%
    dplyr::arrange(desc(get(x)))

  dd <- dplyr::mutate(dd, Description = factor(.data$Description, rev(.data$Description)))

  p <- ggplot(dd) +
    aes(x = .data[[x]], y = .data[['Description']], fill = .data[[color.by]]) +
    geom_bar(stat = "identity", width = bar.width, alpha = 0.7, color = ifelse(add.bar.border, "black", NA)) +
    labs(fill = color.title, x = x.lab, y = NULL, title = title) +
    scale_fill_gradientn(colours = colors)

  if (y.label.position == 'on') {
    p <- p +
      geom_text(aes(x = max(.data[[x]])/100, label = .data[['Description']]), hjust = 0)
  } else {
    p <- p +
      scale_y_discrete(labels = Hmisc::capitalize, position = y.label.position)
  }

  p <- p +
    ggtheme +
    theme(axis.text.x = element_text(size = 10, colour = "black"),
          axis.title.x = element_text(size = 13, colour = "black", face = "bold"),
          axis.text.y = element_text(size = 13, colour = "black"),
          axis.title = element_blank(),
          panel.background = element_blank(),
          legend.text = element_text(size = 11, colour = "black"),
          legend.title = element_text(size = 13, colour = "black", face = "bold"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.position = legend.position,
          plot.title = element_text(hjust = 0.5, size = 14, colour = "black", face = "bold"),
          ...)
  if (y.label.position == 'on') {
    p <- p +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
  }
  return(p)
}



#' KEGG pathway graph visualization
#' @description
#' Simple visualization of KEGG pathway based on pathview package.
#'
#' @param gene.data The same in `pathview::pathview()`.
#' @param gene.type The same in `pathview::pathview()`.
#' @param pathway.id The same in `pathview::pathview()`.
#' @param species The same in `pathview::pathview()`.
#' @param figure.suffix The same in `pathview::pathview()`.
#' @param save.dir The dir path to save the picture.
#'
#' @importFrom pathview pathview
#' @importFrom png readPNG
#' @importFrom graphics par plot.new rasterImage
#'
#' @export
#'
simple_KEGG_graph <- function(gene.data, gene.type = "SYMBOL", pathway.id,
                              species = "hsa", figure.suffix, save.dir = getwd()){
  file_path <- file.path(save.dir, paste0(gsub("\\D", "", pathway.id), ".", figure.suffix, ".png"))
  cli::cli_alert_info('Saving picture to {.file {file_path}} ...')
  pathview::pathview(gene.data = gene.data, gene.idtype = gene.type,
                     pathway.id = gsub("\\D", "", pathway.id),
                     kegg.dir = save.dir, species = species,
                     kegg.native = T, out.suffix = figure.suffix,
                     low = list(gene = "#0a9396", cpd = "blue"),
                     mid = list(gene = "#e9d8a6", cpd = "gray"),
                     high = list(gene = "#bb3e03", cpd = "yellow"))
  imgpng <- png::readPNG(file_path)
  graphics::par(mar = c(0, 0, 0, 0))
  graphics::plot.new()
  graphics::rasterImage(imgpng, 0, 0, 1, nrow(imgpng)/ncol(imgpng))

}


