#' runNMF
#'
#' perform NMF
#'
#' @param object A Seurat object
#' @param savePath A path to save the results files(suggest to create a foler named by sample name)
#' @param rank A integer value to determine the dimension of decomposition matrix
#'
#'
#' @return A list containing a seuratObject, a list of results of intermediate process,
#' a list of plots and a logical value to show whether the process completed successfully
#' @export
#'
#' @import Seurat ggplot2
#' @importFrom NNLM nnmf
#' @importFrom pheatmap pheatmap
#' @importFrom stats median quantile
#' @importFrom cowplot plot_grid
#'
runNMF <- function(object,
                   savePath,
                   rank = 5,
                   base.size = 8,
                   pt.size = 1.6){

  # collect results and plots
  results.collector <- list()
  plots.collector <- list()
  bool.completed <- FALSE

  prog.colors <- colorRampPalette(c("white", "yellow", "red"))(50)

  x.data <- GetAssayData(object, slot = "data")

  nmf.results <- NNLM::nnmf(as.matrix(x.data), k = rank, verbose = F)

  colnames(nmf.results$W) <- paste0("p", 1:rank)
  rownames(nmf.results$H) <- paste0("p", 1:rank)
  colnames(nmf.results$H) <- colnames(x.data)

  saveRDS(nmf.results,
          file = file.path(savePath, paste0("NMF-", rank, ".RDS")))

  p.W <- nmf.results$W
  p.H <- nmf.results$H

  results.collector[["nmf.results.H"]] <- p.H
  results.collector[["nmf.results.W"]] <- p.W

  up.lim <- stats::quantile(p.H, 0.99)
  p.H[p.H > up.lim] <- up.lim

  pheatmap::pheatmap(p.H,
                     show_colnames = F,
                     width = 6,
                     height = 2.5,
                     treeheight_row = 10,
                     treeheight_col = 14,
                     filename = file.path(savePath, paste0("NMF-", rank, "-heatmap.png")))

  plots.collector[[paste0("NMF-", rank, "-heatmap")]] <- pheatmap::pheatmap(p.H,
                                                                            show_colnames = F,
                                                                            width = 6,
                                                                            height = 2.5,
                                                                            treeheight_row = 10,
                                                                            treeheight_col = 14)

  p.H <- t(p.H)
  mat.colnames <- make.names(colnames(p.H))
  for(i in 1 : length(mat.colnames)){
    object <- AddMetaData(object, p.H[, i], col.name = mat.colnames[[i]])
  }

  ps <- list()

  for(p.i in  1 : length(mat.colnames)){
    ps[[p.i]] <- Spatial_Plot(object,
                              feature = mat.colnames[[p.i]],
                              show.tissue = F,
                              crop = T,
                              title = colnames(p.H)[[p.i]],
                              legend.title = "",
                              legend.position = "none",
                              colors = prog.colors,
                              discrete = F,
                              base.size = base.size,
                              pt.size = pt.size)
  }

  plots.collector[["p"]] <- ps

  n.col = getPlotCol(ncol(p.H))
  p.all <- cowplot::plot_grid(plotlist = ps,
                              ncol = n.col)
  ggplot2::ggsave(file.path(savePath, paste0("NMF-", rank, "-sp.png")),
                  p.all,
                  width = n.col*4.2,
                  height = ceiling(length(ps) / n.col) * 4.9,
                  dpi = 300)

  plots.collector[[paste0("NMF-", rank, "-sp")]] <- p.all

  bool.completed <- TRUE
  return(list(object = object,
              results = results.collector,
              plots = plots.collector,
              bool.completed = bool.completed))
}
