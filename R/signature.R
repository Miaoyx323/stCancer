#' cellTypeScore
#'
#' calculate the score of some cell types
#'
#' @param SeuratObject A Seurat object
#' @param geneSets
#' @param select
#' @param assay
#' @param method
#' @param savePath
#'
#'
#' @return A Seurat object after scoring
#' @export
#'
#' @import Seurat ggplot2 GSVA
#'

cellTypeScore <- function(SeuratObject,
                          savePath = NULL,
                          geneSets = NULL,
                          select = "main",
                          assay = "SCT",
                          method = "average",
                          species = "human",
                          total.image.name = "cellTypeScore.png",
                          ...){

  # collect results and plots
  results.collector <- list()
  plots.collector <- list()
  bool.completed <- FALSE

  if(is.null(geneSets)){
    ct.anno <- read.table(file.path(system.file(package = "stCancer"), "txt/ct.anno.txt"),
                          sep = "\t",
                          header = T,
                          stringsAsFactors = F)
    if(select == "main"){
      ct.anno <- ct.anno %>% subset(type == "Main")
    }
    geneSets <- list()
    for(i in 1 : nrow(ct.anno)){
      geneSets[ct.anno$name[i]] <- strsplit(gsub(" ", "", ct.anno$genes[i]), split = ",")
    }
  }

  if(species == "mouse"){
    for(i in 1 : length(geneSets)){
      geneSets[[i]] <- getMouseGene(geneSets[[i]])
    }
  }

  data <- GetAssayData(SeuratObject, assay = assay)
  data <- t(data)

  if(method == "AddModuleScore"){
    t.SeuratObject <- AddModuleScore(SeuratObject,
                                     features = geneSets,
                                     assay = "SCT",
                                     name = "geneSet")
    t.scores <- t.SeuratObject[[paste0("geneSet", 1:length(geneSets))]]
    t.scores <- scale(t.scores)
  }else if(method == "average"){
    t.scores <- data.frame(row.names = rownames(data))
    for(i in 1 : length(geneSets)){
      com.genes <- intersect(unlist(geneSets[i]), colnames(data))
      t.scores[[names(geneSets)[i]]] <- rowMeans(as.data.frame(data[, com.genes]))
    }
    t.scores <- as.matrix(t.scores)
  }else if(method == "GSVA"){
    tmp.data <- as.matrix(GetAssayData(SeuratObject,
                                       slot = "scale.data"))
    tmp.data <- tmp.data[VariableFeatures(SeuratObject), ]
    t.scores <- t(gsva(tmp.data, geneSets))
  }else{
    stop(paste0("CellTypeScore: No method called ", method))
  }

  colnames(t.scores) <- names(geneSets)
  colnames(t.scores) <- sapply(colnames(t.scores), function(x){
    gsub(" ", "\\.", x)
  })

  # t.scores <- t.scores[, colSums(is.na(t.scores)) == 0]
  t.scores[is.na(t.scores)] <- 0

  if(ncol(t.scores) == 0){
    return(list(object = SeuratObject,
                results = results.collector,
                plots = plots.collector,
                bool.completed = bool.completed))
  }

  SeuratObject <- AddMetaData(SeuratObject,
                              as.data.frame(t.scores))

  results.collector[["scores"]] <- t.scores

  cells <- colnames(t.scores)
  plot_list <- list()
  plot_dr_list <- list()

  for(i in 1 : length(cells)){
    cellType <- cells[[i]]
    plot_list[[i]] <- ggplotGrob(Spatial_Plot(SeuratObject,
                                              cellType,
                                              discrete = F,
                                              crop = T,
                                              pt.size = 1.6,
                                              base.size = 8,
                                              legend.title = "score",
                                              title = cellType,
                                              ...))
    if(!is.null(savePath)){
      ggsave(file.path(savePath, paste0(cellType, ".png")),
             plot_list[[i]],
             dpi = 300,
             height = 5,
             width = 6)
    }

    plot_dr_list[[i]] <- ggplotGrob(pointDRPlot(SeuratObject,
                                                cellType,
                                                discrete = F,
                                                legend.title = "score",
                                                title = cellType,
                                                ...))

    if(!is.null(savePath)){
      ggsave(file.path(savePath, paste0(cellType, "_dr.png")),
             plot_dr_list[[i]],
             dpi = 300,
             height = 5,
             width = 6)
    }
  }
  plots.collector[["cell.types"]] <- plot_list
  plots.collector[["cell.types.dr"]] <- plot_dr_list

  total_plot <- plot_grid(plotlist = plot_list, ncol = getPlotCol(length(plot_list)))
  total_dr_plot <- plot_grid(plotlist = plot_dr_list, ncol = getPlotCol(length(plot_dr_list)))

  plots.collector[["total.plot"]] <- total_plot
  plots.collector[["total.plot.dr"]] <- total_dr_plot

  if(!is.null(total.image.name) & !is.null(savePath)){
    ggsave(file.path(savePath, total.image.name),
           total_plot,
           dpi = 300,
           height = ceiling(length(plot_list) / getPlotCol(length(plot_list))) * 4.8,
           width = getPlotCol(length(plot_list)) * 4.8)

    ggsave(file.path(savePath, paste0("dr_", total.image.name)),
           total_dr_plot,
           dpi = 300,
           height = ceiling(length(plot_dr_list) / getPlotCol(length(plot_dr_list))) * 4.5,
           width = getPlotCol(length(plot_dr_list)) * 5)
  }

  bool.completed <- TRUE
  return(list(object = SeuratObject,
              results = results.collector,
              plots = plots.collector,
              bool.completed = bool.completed))
}


#' CellCycleAnalysis
#'
#' calculate the score of cell cycle
#'
#' @param SeuratObject A Seurat object
#' @param select
#' @param assay
#' @param method
#' @param savePath
#'
#'
#' @return A Seurat object after NMF
#' @export
#'
#' @import Seurat ggplot2 GSVA
#'
CellCycleAnalysis <- function(SeuratObject,
                              savePath = NULL,
                              geneSets = NULL,
                              ...){

  # collect results and plots
  results.collector <- list()
  plots.collector <- list()
  bool.completed <- FALSE

  # sampleName <- SeuratObject@project.name

  # if(!is.null(savePath)){
  #   savePath <- R.utils::getAbsolutePath(savePath)
  #
  #   savePath_basic <- file.path(savePath, sampleName)
  #   savePath_pheno <<- file.path(savePath_basic, "phenotype")
  #
  #   if(!dir.exists(savePath_pheno)){
  #     dir.create(savePath_pheno, recursive = T)
  #   }
  # }

  if(is.null(geneSets)){
    cycle.marker <- read.table(file.path(system.file(package = "stCancer"), "txt/cellCycle-genes.txt"),
                               header = F,
                               stringsAsFactors = F)
    geneSets <- list()
    geneSets["CellCycle"] <- as.list(cycle.marker)
  }

  results <- cellTypeScore(SeuratObject = SeuratObject,
                           savePath = NULL,
                           geneSets = geneSets,
                           total.image.name = NULL,
                           colors = c("white", "#009b45"),
                           ...)

  SeuratObject <- results$object
  results.collector <- results$results
  bool.completed <- results$bool.completed

  if(is.null(SeuratObject@meta.data[["CellCycle"]])){
    return(list(object = SeuratObject,
                results = results.collector,
                plots = plots.collector,
                bool.completed = bool.completed))
  }

  p <- Spatial_Plot(SeuratObject,
                    "CellCycle",
                    discrete = F,
                    crop = T,
                    base.size = 8,
                    pt.size = 1.6,
                    colors = c("white", "#009b45"))

  plots.collector[["CellCycle"]] <- p

  if(!is.null(savePath)){
    ggsave(file.path(savePath, "cellCycle.png"),
           p,
           dpi = 300,
           height = 5,
           width = 6)
  }

  p <- pointDRPlot(SeuratObject,
                   "CellCycle",
                   discrete = F,
                   legend.title = "score",
                   title = "CellCycle",
                   colors = c("white", "#009b45"))

  plots.collector[["CellCycle_dr"]] <- p

  if(!is.null(savePath)){
    ggsave(file.path(savePath, "cellCycle_dr.png"),
           p,
           dpi = 300,
           height = 5,
           width = 6)
  }

  return(list(object = SeuratObject,
              results = results.collector,
              plots = plots.collector,
              bool.completed = bool.completed))
}

#' StemnessAnalysis
#'
#' Estimate cell stemness according to the Spearman correlation with stemness signature.
#'
#' @param SeuratObject A seurat object to estimate stemness.
#' @param stem.sig An array of stemness signature. The default is NULL, and a prepared signature will be used.
#'
#' @importFrom stats cor
#'
#' @return
#' @export
#'
StemnessAnalysis <- function(SeuratObject,
                             savePath = NULL,
                             geneSets = NULL,
                             ...){
  # collect results and plots
  results.collector <- list()
  plots.collector <- list()
  bool.completed <- FALSE

  if(is.null(geneSets)){
    stem.marker <- read.table(file.path(system.file(package = "stCancer"), "txt/stem.txt"),
                              header = F,
                              stringsAsFactors = F)
    geneSets <- list()
    geneSets["stem"] <- as.list(stem.marker)
  }

  results <- cellTypeScore(SeuratObject = SeuratObject,
                           savePath = NULL,
                           geneSets = geneSets,
                           total.image.name = NULL,
                           colors = c("white", "#ff9000"),
                           ...)
  SeuratObject <- results$object
  results.collector <- results$results
  bool.completed <- results$bool.completed

  if(is.null(SeuratObject@meta.data[["stem"]])){
    return(list(object = SeuratObject,
                results = results.collector,
                plots = plots.collector,
                bool.completed = bool.completed))
  }

  p <- Spatial_Plot(SeuratObject,
                    "stem",
                    discrete = F,
                    crop = T,
                    pt.size = 1.6,
                    base.size = 8,
                    colors = c("white", "#ff9000"))

  plots.collector[["stem"]] <- p

  if(!is.null(savePath)){
    ggsave(file.path(savePath, "stem.png"),
           p,
           dpi = 300,
           height = 5,
           width = 6)
  }


  p <- pointDRPlot(SeuratObject,
                   "stem",
                   discrete = F,
                   legend.title = "score",
                   title = "stem",
                   colors = c("white", "#ff9000"))

  plots.collector[["stem_dr"]] <- p

  if(!is.null(savePath)){
    ggsave(file.path(savePath, "stem_dr.png"),
           p,
           dpi = 300,
           height = 5,
           width = 6)
  }

  return(list(object = SeuratObject,
              results = results.collector,
              plots = plots.collector,
              bool.completed = bool.completed))
}


#' EMTAnalysis
#'
#' calculate the score of EMT
#'
#' @param SeuratObject A Seurat object
#' @param select
#' @param assay
#' @param method
#' @param savePath
#'
#'
#' @return A Seurat object after NMF
#' @export
#'
#' @import Seurat ggplot2 GSVA
#'
EMTAnalysis <- function(SeuratObject,
                        savePath = NULL,
                        geneSets = NULL,
                        ...){

  # collect results and plots
  results.collector <- list()
  plots.collector <- list()
  bool.completed <- FALSE

  if(is.null(geneSets)){
    EMT.marker <- read.table(file.path(system.file(package = "stCancer"), "txt/EMT.txt"),
                             header = F,
                             stringsAsFactors = F)
    geneSets <- list()
    geneSets["EMT"] <- as.list(EMT.marker)
  }

  results <- cellTypeScore(SeuratObject = SeuratObject,
                           savePath = NULL,
                           geneSets = geneSets,
                           total.image.name = NULL,
                           colors = c("white", "#990099"),
                           ...)
  SeuratObject <- results$object
  results.collector <- results$results
  bool.completed <- results$bool.completed

  if(is.null(SeuratObject@meta.data[["EMT"]])){
    return(list(object = SeuratObject,
                results = results.collector,
                plots = plots.collector,
                bool.completed = bool.completed))
  }


  p <- Spatial_Plot(SeuratObject,
                    "EMT",
                    discrete = F,
                    crop = T,
                    pt.size = 1.6,
                    base.size = 8,
                    colors = c("white", "#990099"))

  plots.collector[["EMT"]] <- p

  if(!is.null(savePath)){
    ggsave(file.path(savePath_pheno, "EMT.png"),
           p,
           dpi = 300,
           height = 5,
           width = 6)
  }

  p <- pointDRPlot(SeuratObject,
                   "EMT",
                   discrete = F,
                   legend.title = "score",
                   title = "EMT",
                   colors = c("white", "#990099"))

  plots.collector[["EMT_dr"]] <- p

  if(!is.null(savePath)){
    ggsave(file.path(savePath_pheno, "EMT_dr.png"),
           p,
           dpi = 300,
           height = 5,
           width = 6)
  }

  return(list(object = SeuratObject,
              results = results.collector,
              plots = plots.collector,
              bool.completed = bool.completed))
}


#' TLSAnalysis
#'
#' calculate the score of TLS
#'
#' @param SeuratObject A Seurat object
#' @param select
#' @param assay
#' @param method
#' @param savePath
#'
#'
#' @return A Seurat object after NMF
#' @export
#'
#' @import Seurat ggplot2 GSVA
#'
TLSAnalysis <- function(SeuratObject,
                        savePath = NULL,
                        geneSets = NULL,
                        ...){

  # collect results and plots
  results.collector <- list()
  plots.collector <- list()
  bool.completed <- FALSE

  if(is.null(geneSets)){
    tls.marker <- read.table(file.path(system.file(package = "stCancer"), "txt/tls.txt"),
                             header = F,
                             stringsAsFactors = F)
    geneSets <- list()
    geneSets["TLS"] <- as.list(tls.marker)
  }

  results <- cellTypeScore(SeuratObject = SeuratObject,
                           savePath = NULL,
                           geneSets = geneSets,
                           total.image.name = NULL,
                           ...)
  SeuratObject <- results$object
  results.collector <- results$results
  bool.completed <- results$bool.completed

  if(is.null(SeuratObject@meta.data[["TLS"]])){
    return(list(object = SeuratObject,
                results = results.collector,
                plots = plots.collector,
                bool.completed = bool.completed))
  }

  p <- Spatial_Plot(SeuratObject,
                    "TLS",
                    discrete = F,
                    crop = T,
                    pt.size = 1.6,
                    base.size = 8,
                    colors = c("white", "#660099"))

  plots.collector[["TLS"]] <- p

  if(!is.null(savePath)){
    ggsave(file.path(savePath_Score, "TLS.png"),
           p,
           dpi = 300,
           height = 5,
           width = 6)
  }

  p <- pointDRPlot(SeuratObject,
                   "TLS",
                   discrete = F,
                   legend.title = "score",
                   title = "TLS",
                   colors = c("white", "#660099"))

  plots.collector[["TLS_dr"]] <- p

  if(!is.null(savePath)){
    ggsave(file.path(savePath_Score, "TLS_dr.png"),
           p,
           dpi = 300,
           height = 5,
           width = 6)
  }

  return(list(object = SeuratObject,
              results = results.collector,
              plots = plots.collector,
              bool.completed = bool.completed))

}
