#' geom_spatial
#'
#' Draw spatial figure
#'
#' @return
#'
#' @import ggplot2
#'
geom_spatial <-  function(mapping = NULL,
                          data = NULL,
                          stat = "identity",
                          position = "identity",
                          na.rm = FALSE,
                          show.legend = NA,
                          inherit.aes = FALSE,
                          ...) {
  GeomCustom <- ggproto(
    "GeomCustom",
    Geom,
    setup_data = function(self, data, params) {
      data <- ggproto_parent(Geom, self)$setup_data(data, params)
      data
    },
    draw_group = function(data, panel_scales, coord) {
      vp <- grid::viewport(x=data$x, y=data$y)
      g <- grid::editGrob(data$grob[[1]], vp=vp)
      ggplot2:::ggname("geom_spatial", g)
    },
    required_aes = c("grob","x","y")
  )

  layer(
    geom = GeomCustom,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}


#' histPlot
#'
#' Draw histogram of Seurat object
#'
#' @return
#'
#'
histPlot <- function(object,
                     feature,
                     color = "#a788ab",
                     xlines = c(),
                     y.label = "Spot number"){
  return(histPlot.data.frame(data = object@meta.data,
                             value = feature,
                             color = color,
                             xlines = xlines,
                             y.label = y.label))
}

#' histPlot.data.frame
#'
#' Draw histogram of data.frame
#'
#' @return
#'
#'
histPlot.data.frame <- function(data,
                                value,
                                color = "#a788ab",
                                xlines = c(),
                                y.label = "Spot number"){
  p <- ggplot(data, aes(x = .data[[value]])) +
    geom_histogram(bins = 200, position = "stack", fill = color, alpha = 0.6) +
    labs(x = value, y = y.label) +
    ggplot_config(base.size = 6)

  for(x in xlines){
    p <- p + geom_vline(xintercept = x, colour = "red", linetype = "dashed")
  }
  return(p)
}

#' ggplot_config
#'
#' Default ggplot config
#'
#' @return
#'
#' @import ggplot2
#'
ggplot_config <- function(base.size = 8){
  p <- theme_classic() +
    theme(plot.title = element_text(size = 2 * base.size),
          axis.title.x = element_text(size = 2 * base.size, vjust = -0.2),
          axis.title.y = element_text(size = 2 * base.size, vjust = 0.2),
          axis.text.x = element_text(size = 1.6 * base.size),
          axis.text.y = element_text(size = 1.6 * base.size),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_text(size = 2 * base.size - 1),
          legend.text = element_text(size = 2 * base.size - 1)
    )
  return(p)
}


#' Spatial_Plot
#'
#' The most fundamental function
#'
#' @param object Seurat object
#' @param feature Feature in Seurat meta.data
#'
#' @return
#' @export
#'
#' @import ggplot2
#' @importFrom grid rasterGrob
#' @importFrom dplyr tibble
#' @importFrom dplyr %>%
#'
Spatial_Plot <- function(object,
                         feature,
                         image.key = "slice1",
                         show.tissue = F,
                         crop = F,
                         title = NULL,
                         legend.title = NULL,
                         legend.position = "right",
                         legend.size = 3,
                         colors = NULL,
                         limit.quantile = 0,
                         limit.type = "both",
                         discrete = F,
                         margin = 5,
                         base.size = 8,
                         pt.size = "auto",
                         alpha = 1.0){

  if(!is.null(feature)){
    value <- object@meta.data[[feature]]
  }else{
    value <- NULL
  }
  tissue.positions <- object@images[[image.key]]@coordinates
  scale.factors <- object@images[[image.key]]@scale.factors
  im <- object@images[[image.key]]@image

  if(crop){
    imrow.min <- floor(scale.factors$crop.row.min * scale.factors$lowres) - margin
    imrow.max <- ceiling(scale.factors$crop.row.max * scale.factors$lowres) + margin
    imcol.min <- floor(scale.factors$crop.col.min * scale.factors$lowres) - margin
    imcol.max <- ceiling(scale.factors$crop.col.max * scale.factors$lowres) + margin
    tissue.positions$imagerow <- tissue.positions$imagerow * scale.factors$lowres - imrow.min
    tissue.positions$imagecol <- tissue.positions$imagecol * scale.factors$lowres - imcol.min
    im <- im[imrow.min:imrow.max, imcol.min:imcol.max, ]
    im.grob <- rasterGrob(im, width = unit(1, "npc"), height = unit(1,"npc"))
    im.tibble <- tibble(grob = list(im.grob))
  }else{
    tissue.positions$imagerow <- tissue.positions$imagerow * scale.factors$lowres
    tissue.positions$imagecol <- tissue.positions$imagecol * scale.factors$lowres
    im.grob <- rasterGrob(im, width = unit(1, "npc"), height = unit(1,"npc"))
    im.tibble <- tibble(grob = list(im.grob))
  }

  if(is.null(legend.title)){
    legend.title <- feature
  }

  if(is.null(colors) & !is.null(feature)){
    if(discrete){
      colors <- getDefaultColors(length(unique(value)))
    }else{
      colors <- c("white", "red")
    }
  }else if(all(colors == "cluster.color")){
    tmp <- object@meta.data[["cluster.color"]]
    if(!is.null(tmp)){
      # names(tmp) <- as.numeric(as.character(object@meta.data[["seurat_clusters"]]))
      colors <- unique(tmp)
      names(colors) <- unique(as.numeric(as.character(object@meta.data[["seurat_clusters"]])))
    }else{
      if(discrete){
        colors <- getDefaultColors(length(unique(value)))
      }else{
        colors <- c("white", "red")
      }
    }
  }

  if(!is.null(feature)){
    if(!discrete){
      fill.value <- as.numeric(value)
      dn.q <- ifelse(limit.type %in% c("down", "both"), limit.quantile, 0)
      up.q <- ifelse(limit.type %in% c("up", "both"), 1 - limit.quantile, 1)
      pos.values <- value[value > 0]
      if(length(unique(is.na(pos.values))) == 2){
        dn.thres <- quantile(pos.values, dn.q)
        up.thres <- quantile(pos.values, up.q)
        fill.value <- ifelse(fill.value < dn.thres, dn.thres,
                             ifelse(fill.value > up.thres, up.thres, fill.value))
        # value <- factor(fill.value)
        value <- fill.value
      }else{

      }
    }else{
      value <- factor(value)
    }
  }

  # plot
  if(!is.null(feature)){
    tissue.positions[[feature]] <- value
  }
  p <- tissue.positions %>% ggplot()

  if(show.tissue){
    p <- p + geom_spatial(data = im.tibble, aes(grob = grob), x = 0.5, y = 0.5)
  }

  y.max <- nrow(im)
  x.max <- ncol(im)

  if(pt.size == "auto"){
    if(crop){
      pt.size <- 1.6
    }else{
      pt.size <- 1.2
    }
  }

  if(!is.null(feature)){
    p <- p + geom_point(aes(x = imagecol, y = y.max - imagerow, color = .data[[feature]]),
                        shape = 16,
                        size = pt.size,
                        alpha = alpha)
    # p <- p + geom_point(aes(x = imagecol, y = y.max - imagerow, color = .data[[feature]], alpha = .data[[feature]]),
    #                     shape = 16,
    #                     size = pt.size)
  }

  # setting
  p <- p +
    scale_x_continuous(expand = c(0, 0), limits = c(0, ncol(im))) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, nrow(im))) +
    coord_fixed(ratio = 1) +
    xlab("") + ylab("") +
    ggtitle(title) +
    labs(fill = legend.title) +
    theme_test() + # boundary of the image
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size = 4 * base.size + 2, hjust = 0),
          legend.title = element_text(size = 2 * base.size - 1),
          legend.position = legend.position,
          legend.text = element_text(size = 2 * base.size - 1))


  if(!is.null(value)){
    if(discrete){
      p <- p +
        scale_color_manual(name = legend.title,
                           values = colors) +
        guides(fill = guide_legend(override.aes = list(size = legend.size),
                                   keywidth = 0.1,
                                   keyheight = 0.15,
                                   default.unit = "inch"))
    }else{
      p <- p + scale_color_gradientn(name = legend.title,
                                     colours = colors)
    }
  }
  return(p)
}


#' getPointSize
#'
#' Automatic calculate the points' size
#'
#' @return
#'
#'
getPointSize <- function(spm,
                         crop){
  size_scale <- 150

  if(!crop){
    return((5.8 / spm$width[1]) * size_scale)
  }else{
    return((5.8 / spm$width.cropped[1]) * size_scale - 0.3)
  }
}


#' genePropPlot
#'
#' Drow gene proportion
#'
#' @param gene.manifest
#' @param expr.frac
#'
#' @return
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom gridExtra arrangeGrob
#'
genePropPlot <- function(gene.manifest,
                         expr.frac){
  show.num = 100
  if("bg.percent" %in% colnames(gene.manifest)){
    gene.manifest <- gene.manifest[order(gene.manifest$bg.percent, decreasing = T), ]
  }else{
    gene.manifest <- gene.manifest[order(gene.manifest$prop.median, decreasing = T), ]
  }
  gene.show <- head(gene.manifest$Symbol, show.num)

  rownames(gene.manifest) <- gene.manifest$Symbol
  gene.manifest <- gene.manifest[gene.show, ]

  rate.df.plot <- melt(as.data.frame(as.matrix(t(expr.frac[gene.show, ]))), id.vars = NULL)
  rate.df.plot$Annotation <- factor(gene.manifest[as.character(rate.df.plot$variable), ]$Annotation,
                                    levels = c("mitochondrial", "ribosome", "dissociation", "other"), ordered = T)
  rate.df.plot$variable <- factor(rate.df.plot$variable,
                                  levels = rev(gene.show), ordered = TRUE)

  sub.ix <- c(rep("1-50", 50), rep("51-100", 50))
  names(sub.ix) <- gene.show
  rate.df.plot$subPlot <- factor(sub.ix[as.character(rate.df.plot$variable)],
                                 levels = c("1-50", "51-100"), ordered = T)

  if("bg.percent" %in% colnames(gene.manifest)){
    bg.df <- data.frame(ix = c(1:50, 1:50) + 0.1,
                        frac = rev(gene.manifest[gene.show, ]$bg.percent),
                        type = "Background proportion")
    bg.df$subPlot <- factor(sub.ix[as.character(rev(gene.show))],
                            levels = c("1-50", "51-100"), ordered = T)
  }else{
    bg.df <- NULL
  }

  gene.colors <- c(other = "#d56f6d", ribosome = "#f29721", dissociation = "#65a55d", mitochondrial = "#3778bf")
  p <- ggplot() + coord_flip() +
    scale_color_manual(values = "red", name = "") +
    scale_fill_manual(values = gene.colors, name = "Gene:") +
    scale_y_continuous(breaks = c(0.0001, 0.001, 0.01, 0.1, 1),
                       labels = c("0.0001", "0.001", "0.01", "0.1", "1"),
                       trans = 'log10') +
    theme_bw() +
    theme(axis.title.y = element_blank(),
          legend.position = "top",
          legend.title = element_text(face="bold"))

  p1 <- p + geom_boxplot(data = subset(rate.df.plot, subPlot == "1-50"),
                         mapping = aes(x = variable, y = value, fill = Annotation),
                         outlier.size = 0.1, alpha = 0.8) +
    xlab("") + ylab(paste0("Gene proportion of total UMI (1-50)"))

  p2 <- p + geom_boxplot(data = subset(rate.df.plot, subPlot == "51-100"),
                         mapping = aes(x = variable, y = value, fill = Annotation),
                         outlier.size = 0.1, alpha = 0.8) +
    xlab("") + ylab(paste0("Gene proportion of total UMI (51-100)"))

  all.p <- p + geom_boxplot(data = rate.df.plot,
                            mapping = aes(x = variable, y = value, fill = Annotation),
                            outlier.size = 0.1, alpha = 0.8)

  if(!is.null(bg.df)){
    p1 <- p1 + geom_point(data = subset(bg.df, subPlot == "1-50"),
                          aes(x = ix, y = frac, color = type), shape = "*", size = 5)
    p2 <- p2 + geom_point(data = subset(bg.df, subPlot == "51-100"),
                          aes(x = ix, y = frac, color = type), shape = "*", size = 5)
    all.p <- all.p + geom_point(data = bg.df, aes(x = ix, y = frac, color = type), shape = "*", size = 5)
  }

  p <- grid_arrange_shared_legend(p1, p2, all.p = all.p, ncol = 2, nrow = 1)

  return(p)
}



#' grid_arrange_shared_legend
#'
#' grid arrange
#'
#' @param
#'
#' @return
#'
#' @importFrom gridExtra arrangeGrob
#' @importFrom grid unit.c grid.newpage grid.draw
#'
grid_arrange_shared_legend <- function(...,
                                       all.p,
                                       ncol = length(list(...)),
                                       nrow = 1,
                                       position = c("bottom", "right")) {

  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(all.p + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)

  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)

  # return gtable invisibly
  invisible(combined)
}



#' pointDRPlot
#'
#' plot features in a reduced dimensionality space
#'
#' @param object Seurat object
#' @param feature
#'
#' @return
#' @export
#'
#' @importFrom gridExtra arrangeGrob
#'
pointDRPlot <- function(object,
                        feature,
                        reduction.func = "tsne",
                        colors = NULL,
                        title = NULL,
                        discrete = T,
                        limit.quantile = 0,
                        limit.q.nz = FALSE,
                        point.type = 1,
                        point.size = NULL,
                        base.size = 6,
                        legend.position = "right",
                        legend.title = NULL){

  if(is.null(object@reductions[[reduction.func]])){
    stop("Error in 'pointDRPlot': 'reduction.func' ", reduction.func, " is not in Seurat object.\n")
  }
  if(!(point.type %in% c(1, 2))){
    stop("Error in 'pointDRPlot': 'point.type' ", point.type, " is not allowed.\n")
  }

  if(is.null(legend.title)){
    legend.title <- feature
  }
  if(is.null(colors)){
    if(discrete){
      colors <- getDefaultColors(length(unique(object@meta.data[[feature]])))
    }else{
      colors <- c("white", "red")
    }
  }

  coor <- as.data.frame(object@reductions[[reduction.func]]@cell.embeddings)
  feature <- as.data.frame(object@meta.data[[feature]])
  colnames(feature) <- "feature"

  data <- cbind.data.frame(coor, feature)

  ratio <- diff(range(data[, 1])) / diff(range(data[, 2]))

  fill.feature <- data$feature
  if(!discrete){
    if(limit.q.nz){
      cur.features <- data$feature[data$feature > 0]
      low.thres <- max(0, min(cur.features))
      up.thres <- quantile(cur.features, 1 - limit.quantile)
    }else{
      low.thres <- quantile(data$feature, limit.quantile)
      up.thres <- quantile(data$feature, 1 - limit.quantile)
    }
    fill.feature <- ifelse(fill.feature < low.thres, low.thres,
                           ifelse(fill.feature > up.thres, up.thres, fill.feature))
  }

  p <- ggplot()

  if(point.type == 1){
    if(is.null(point.size)){ point.size <- 1 }
    p <- p +
      geom_point(data,
                 mapping = aes(x = data[, 1],
                               y = data[, 2],
                               fill = fill.feature),
                 shape = 21, size = point.size, stroke = 0.2, color = "lightgrey") +
      coord_fixed(ratio = ratio) +
      ggtitle(title) +
      ggplot_config(base.size = base.size) +
      labs(x = colnames(coor)[1], y = colnames(coor)[2]) +
      labs(fill = legend.title) +
      theme(legend.position = legend.position,
            plot.title = element_text(size = 4 * base.size + 2, hjust = 0),)

    if(discrete){
      p <- p + scale_fill_manual(values = colors,
                                 guide = guide_legend(override.aes = list(size = 3),
                                                      keywidth = 0.1,
                                                      keyheight = 0.15,
                                                      default.unit = "inch"))
    }else{
      p <- p + scale_fill_gradientn(colors = colors)
    }
  }else if(point.type == 2){
    if(is.null(point.size)){ point.size <- 0.2 }
    p <- p + geom_point(data[, ],
                        mapping = aes(x = data[, 1],
                                      y = data[, 2],
                                      color = fill.feature),
                        shape = 16, size = point.size) +
      coord_fixed(ratio = ratio) +
      ggtitle(title) +
      ggplot_config(base.size = base.size) +
      labs(x = colnames(coor)[1], y = colnames(coor)[2]) +
      labs(color = legend.title) +
      theme(legend.position = legend.position,
            plot.title = element_text(size = 4 * base.size + 2, hjust = 0),)

    if(discrete){
      p <- p + scale_color_manual(values = colors,
                                  guide = guide_legend(override.aes = list(size = 3),
                                                       keywidth = 0.1,
                                                       keyheight = 0.15,
                                                       default.unit = "inch"))
    }else{
      p <- p + scale_color_gradientn(colors = colors)
    }
  }
  return(p)
}


#' preDEheatmap
#'
#' Do preparation for heatmap
#'
#' @param object Seurat object
#'
#' @return
#' @export
#'
#' @import Seurat
#'
preDEheatmap <- function(object,
                         col.name = "seurat_clusters",
                         genes = NULL,
                         spots = NULL,
                         assay = "SCT",
                         slot = "scale.data",
                         min.value = -2.5,
                         max.value = 2.5){

  expr.data <- GetAssayData(object = object, slot = slot, assay = assay)
  if(!is.null(genes)){
    genes <- intersect(genes, rownames(expr.data))
    expr.data <- expr.data[genes, ]
  }
  if(!is.null(spots)){
    spots <- intersect(spots, colnames(expr.data))
    expr.data <- expr.data[, spots]
  }

  spot.cluster <- as.data.frame(object@meta.data[[col.name]])
  # spot.cluster <- as.data.frame(paste0("cluster", as.character(object@meta.data[[col.name]])))
  rownames(spot.cluster) <- colnames(object)
  colnames(spot.cluster) <- col.name

  ind <- order(spot.cluster[[col.name]])
  spot.cluster <- as.data.frame(spot.cluster[ind, ])
  rownames(spot.cluster) <- colnames(object)[ind]
  colnames(spot.cluster) <- col.name

  expr.data <- expr.data[, rownames(spot.cluster)]

  ## limitData
  expr.data <- limitData(expr.data, min = min.value, max = max.value)

  ## gaps_col
  num.cluster <- table(spot.cluster[[col.name]])
  num.cluster <- num.cluster[as.character(min(as.numeric(names(num.cluster))) :
                                            max(as.numeric(names(num.cluster))))]
  gaps_col <- cumsum(num.cluster)

  return(list(expr.data = expr.data,
              spot.cluster = spot.cluster,
              gaps_col = gaps_col))
}



#' getPlotCol
#'
#' Get plots col.num
#'
#' @param num
#'
#' @return
#'
getPlotCol <- function(num){
  num_sqrt <- floor(sqrt(num))
  col.num <- ceiling(num / num_sqrt)
  return(ifelse(col.num >= 7, 7, col.num))
}


#' clusterBarPlot
#'
#' Draw
#'
#' @param spot.manifest
#'
#' @return
#'
clusterBarPlot <- function(spot.annotation,
                           sel.col = "spot.Type",
                           spot.colors = NULL,
                           legend.title = NULL,
                           legend.position = "bottom"){
  if(is.null(legend.title)){
    legend.title <- sel.col
  }
  if(is.null(spot.colors)){
    spot.colors <- getDefaultColors(length(unique(spot.annotation[[sel.col]])))
  }

  bar.df <- melt(table(spot.annotation[c(sel.col, "Cluster")]))
  # bar.df$Cluster = factor(bar.df$Cluster,
  #                         levels = 0:(length(unique(bar.df$Cluster))-1))
  bar.df$Cluster = factor(bar.df$Cluster)

  p <- ggplot(bar.df, aes(x = Cluster, y = value, fill = .data[[sel.col]])) +
    geom_bar(stat = "identity") +
    ggplot_config(base.size = 6) +
    scale_fill_manual(values = spot.colors) +
    labs(y = "Number of spots") +
    guides(fill = guide_legend(override.aes = list(size=1), title = legend.title)) +
    theme(legend.position = legend.position)

  return(p)
}


#' GeneSpatialPlot
#'
#' Plot gene expression level in spatial
#'
#' @param object Seurat object
#' @param gene gene in object
#'
#' @return
#' @export
#'
#' @import ggplot2
#' @importFrom grid rasterGrob
#' @importFrom dplyr tibble
#' @importFrom dplyr %>%
#'
GeneSpatialPlot <- function(object,
                            gene,
                            ...){
  tol.genes <- rownames(object)
  if(!(gene %in% tol.genes)){
    stop("No gene called ", gene)
  }

  object <- AddMetaData(object,
                        as.data.frame(t(as.matrix(GetAssayData(object[gene, ])))),
                        col.name = make.names(gene))

  return(Spatial_Plot(object = object,
                      feature = make.names(gene),
                      discrete = F,
                      ...))
}
