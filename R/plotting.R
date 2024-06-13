SpatialFeature_Plot_single <- function(object,
                                       feature,
                                       crop = TRUE,
                                       stroke = NA,
                                       legend.name = NULL,
                                       legend.color = getDefaultFeatureColors(type = 'seq'),
                                       pt.size.factor = 1.55,
                                       ...){
    suppressMessages(suppressWarnings(
        p <- SpatialFeaturePlot(object,
                                features = feature,
                                stroke = NA,
                                crop = crop,
                                pt.size.factor = pt.size.factor,
                                ...) +
            theme(legend.position = 'right') +
            scale_fill_gradientn(
                name = ifelse(is.null(legend.name),
                              feature,
                              legend.name),
                colours = legend.color(100)
            )
    ))
    return(p)
}

#' SpatialFeature_Plot
#'
#' Draw continuous features in the spatial location
#'
#' @return plots
#'
#' @import ggplot2
#' @import patchwork
#'
SpatialFeature_Plot <- function(object,
                                features,
                                crop = TRUE,
                                stroke = NA,
                                legend.name = NULL,
                                legend.color = getDefaultFeatureColors(type = 'seq'),
                                ...) {
    if(length(features) == 1){
        return(SpatialFeature_Plot_single(object,
                                          features,
                                          crop = crop,
                                          stroke = stroke,
                                          legend.name = legend.name,
                                          legend.color = legend.color,
                                          ...))
    }else{
        plots <- vector(mode = 'list', length = length(features))
        for(i in 1:length(features)){
            plots[[i]] <- SpatialFeature_Plot_single(object,
                                                     features,
                                                     crop = crop,
                                                     stroke = stroke,
                                                     legend.name = ifelse(is.null(legend.name),
                                                                          features[i],
                                                                          legend.name[i]),
                                                     legend.color = legend.color,
                                                     ...)
        }
        return(wrap_plots(plots = plots))
    }
}


SpatialDim_Plot_single <- function(object,
                                   group.by = NULL,
                                   crop = TRUE,
                                   stroke = NA,
                                   legend.name = NULL,
                                   legend.color = NULL,
                                   pt.size.factor = 1.65,
                                   ...){
    suppressMessages(suppressWarnings(
        p <- SpatialDimPlot(object,
                            group.by = group.by,
                            crop = crop,
                            stroke = stroke,
                            pt.size.factor = pt.size.factor,
                            ...) +
            theme(legend.position = 'right',
                  legend.key = element_blank()) +
            guides(fill=guide_legend(override.aes = list(size=4)))
    ))
    if(!is.null(legend.color)){
        if(!is.null(group.by) | !is.null(legend.name)){
            suppressMessages(suppressWarnings(
                p <- p + scale_fill_manual(name = ifelse(is.null(legend.name),
                                                         group.by,
                                                         legend.name),
                                           values = legend.color)
            ))
        }else{
            suppressMessages(suppressWarnings(
                p <- p + scale_fill_manual(values = legend.color)
            ))
        }
    }
    return(p)
}

#' SpatialDim_Plot
#'
#' Draw discrete features in the spatial location
#'
#' @return plots
#'
#' @import ggplot2
#' @import patchwork
#'
SpatialDim_Plot <- function(object,
                            group.by = NULL,
                            crop = TRUE,
                            stroke = NA,
                            legend.name = NULL,
                            legend.color = NULL,
                            ...) {
    if(is.null(group.by)){
        return(SpatialDim_Plot_single(object,
                                      group.by = group.by,
                                      crop = crop,
                                      stroke = stroke,
                                      legend.name = legend.name,
                                      legend.color = legend.color,
                                      ...))
    }
    if(length(group.by) == 1){
        return(SpatialDim_Plot_single(object,
                                      group.by = group.by,
                                      crop = crop,
                                      stroke = stroke,
                                      legend.name = legend.name,
                                      legend.color = legend.color,
                                      ...))
    }else{
        plots <- vector(mode = 'list', length = length(group.by))
        for(i in 1:length(group.by)){
            plots[[i]] <- SpatialDim_Plot_single(object,
                                                 group.by[i],
                                                 crop = crop,
                                                 stroke = stroke,
                                                 legend.name = ifelse(is.null(legend.name),
                                                                      NULL,
                                                                      legend.name[i]),
                                                 legend.color = ifelse(is.null(legend.color),
                                                                       NULL,
                                                                       legend.color),
                                                 ...)
        }
        return(wrap_plots(plots = plots))
    }
}

Dim_Plot_single <- function(object,
                            group.by = NULL,
                            legend.name = NULL,
                            legend.color = NULL,
                            ...){
    p <- DimPlot(object,
                 group.by = group.by,
                 ...)
    if(is.null(legend.color)){
        return(p)
    }else{
        if(is.null(legend.name)){
            p <- p +
                scale_color_manual(values = legend.color)
        }else{
            p <- p +
                scale_color_manual(name = legend.name,
                                   values = legend.color)
        }
        return(p)
    }
}

Dim_Plot <- function(object,
                     group.by = NULL,
                     legend.name = NULL,
                     legend.color = NULL,
                     ...){
    if(is.null(legend.color)){
        return(DimPlot(object = object,
                       group.by = group.by))
    }

    if(length(group.by) == 1){
        return(Dim_Plot_single(object,
                               group.by = group.by,
                               legend.name = legend.name,
                               legend.color = legend.color,
                               ...))
    }else{
        plots <- vector(mode = 'list', length = length(group.by))
        for(i in 1:length(group.by)){
            plots[[i]] <- Dim_Plot_single(object,
                                          group.by[i],
                                          legend.name = ifelse(is.null(legend.name),
                                                               NULL,
                                                               legend.name[i]),
                                          legend.color = legend.color[i],
                                          ...)
        }
        return(wrap_plots(plots = plots))
    }
}

Feature_Plot_single <- function(object,
                                feature,
                                legend.name = NULL,
                                legend.color = getDefaultFeatureColors(type = 'seq'),
                                ...){
    suppressMessages(suppressWarnings(
        p <- FeaturePlot(object,
                         features = feature,
                         ...) +
            theme(legend.position = 'right') +
            scale_color_gradientn(
                name = ifelse(is.null(legend.name),
                              feature,
                              legend.name),
                colours = legend.color(100)
            )
    ))
    return(p)
}

#' Feature_Plot
#'
#' Draw continuous features
#'
#' @return plots
#'
#' @import ggplot2 patchwork
#'
Feature_Plot <- function(object,
                         features,
                         legend.name = NULL,
                         legend.color = getDefaultFeatureColors(type = 'seq'),
                         ...) {
    if(length(features) == 1){
        return(Feature_Plot_single(object,
                                   features,
                                   legend.name = legend.name,
                                   legend.color = legend.color,
                                   ...))
    }else{
        plots <- vector(mode = 'list', length = length(features))
        for(i in 1:length(features)){
            plots[[i]] <- Feature_Plot_single(object,
                                              features[i],
                                              legend.name = ifelse(is.null(legend.name),
                                                                   features[i],
                                                                   legend.name[i]),
                                              legend.color = legend.color,
                                              ...)
        }
        return(wrap_plots(plots = plots))
    }
}


#' histPlot
#'
#' Draw histogram of Seurat object
#'
#' @return plots
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
#' @return plots
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
#' @return plots
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





#' getPointSize
#'
#' Automatic calculate the points' size
#'
#' @return number
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
#' @param gene.manifest A dataframe of gene information
#' @param expr.frac A matrix including the expression fraction of each gene
#'
#' @return plots
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
    gene.show <- head(rownames(gene.manifest), show.num)

    # rownames(gene.manifest) <- gene.manifest$Symbol
    gene.manifest <- gene.manifest[gene.show, ]

    rate.df.plot <- reshape2::melt(as.data.frame(as.matrix(t(expr.frac[gene.show, ]))), id.vars = NULL)
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
#'
#' @return plots
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
                       "bottom" = gridExtra::arrangeGrob(do.call(gridExtra::arrangeGrob, gl),
                                              legend,
                                              ncol = 1,
                                              heights = grid::unit.c(unit(1, "npc") - lheight, lheight)),
                       "right" = gridExtra::arrangeGrob(do.call(gridExtra::arrangeGrob, gl),
                                             legend,
                                             ncol = 2,
                                             widths = grid::unit.c(unit(1, "npc") - lwidth, lwidth)))
    grid::grid.newpage()
    grid::grid.draw(combined)

    # return gtable invisibly
    invisible(combined)
}



#' preDEheatmap
#'
#' Do preparation for heatmap
#'
#' @param object Seurat object
#'
#' @return plots
#' @export
#'
#' @import Seurat
#'
preDEheatmap <- function(object,
                         col.name = "seurat_clusters",
                         genes = NULL,
                         spots = NULL,
                         assay = NULL,
                         slot = "scale.data",
                         min.value = -2.5,
                         max.value = 2.5){

    if(is.null(assay)){
        assay <- object@active.assay
    }

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
#' @return Number of columns of plots
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
#' @return plots
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

    bar.df <- reshape2::melt(table(spot.annotation[c(sel.col, "Cluster")]))
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

