#' stStatistics
#'
#' perform spot QC, gene QC, visualization.
#'
#' @param sampleName A character string giving a label for this sample.
#' @param dataPath A path containing the cell ranger processed data.
#' Under this path, folders 'filtered_feature_bc_matrix' and 'raw_feature_bc_matrix' exist generally.
#' @param savePath A path to save the results files(suggest to create a foler named by sample name).
#' @param crop A logical value indicating whether to crop the image during plotting
#' @param rm.isolated A logical value indicating whether to remove isolated spots
#' @param region.threshold A integer value, if the area of connected domain of spots are less than `region.threshold`,
#' they will be removed.
#' @param authorName A character string for authors name and will be shown in the report.
#' @param genReport A logical value indicating whether to generate a .html/.md report (suggest to set TRUE).
#'
#'
#' @return A list containing a seuratObject, a list of results of intermediate process and
#' a list of plots
#' @export
#'
#' @import Seurat Matrix ggplot2 knitr
#' @importFrom markdown markdownToHTML
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats median quantile
#' @importFrom utils read.delim read.table write.csv write.table read.csv
#' @importFrom cowplot plot_grid
#'
stStatistics <- function(sampleName,
                         dataPath,
                         savePath,
                         species = "human",
                         h5 = FALSE,
                         crop = TRUE,
                         rm.isolated = TRUE,
                         region.threshold = 3,
                         authorName = NULL,
                         genReport = TRUE,
                         ...){

  message("[", Sys.time(), "] START: RUN stStatistics")

  print(sampleName)

  dataPath <- filePathCheck(dataPath)
  savePath <- filePathCheck(savePath)

  dataPath <- R.utils::getAbsolutePath(dataPath)
  savePath <- R.utils::getAbsolutePath(savePath)

  savePath_basic <- file.path(savePath, sampleName)
  savePath_data <- file.path(savePath_basic, "data/")
  savePath_qc <- file.path(savePath_basic, "QC/")

  if(!dir.exists(savePath)){
    dir.create(savePath, recursive = T)
  }

  if(!dir.exists(savePath_basic)){
    dir.create(savePath_basic, recursive = T)
  }

  if(!dir.exists(savePath_data)){
    dir.create(savePath_data, recursive = T)
  }

  if(!dir.exists(savePath_qc)){
    dir.create(savePath_qc, recursive = T)
  }

  # collect results and plots
  results.collector <- list()
  plots.collector <- list()

  ## ------ Read data ------
  spatial <- Load_Spatial(dataPath,
                          rm.isolated = rm.isolated,
                          region.threshold = region.threshold,
                          sampleName = sampleName,
                          h5 = h5,
                          ...)

  ## ------ Quality control of spot ------
  # remove the spot that contains no barcode
  valid_spot <- colnames(spatial)[which(spatial@meta.data[["nCount_Spatial"]] != 0)]
  spatial <- spatial[, valid_spot]

  # calculate percents of mito.genes and ribo.genes in each spot
  if(species == "human"){
    mito.genes <- grep('^MT-', rownames(spatial), value = TRUE)
    ribo.genes <- grep('^RPL|^RPS|^MRPL|^MRPS', rownames(spatial), value = TRUE)
    ig.genes <- grep('^IG', rownames(spatial), value = TRUE)
  }else if(species == "mouse"){
    mito.genes <- grep('^mt-', rownames(spatial), value = TRUE)
    ribo.genes <- grep('^Rpl|^Rps|^Mrpl|^Mrps', rownames(spatial), value = TRUE)
    ig.genes <- grep('^Ig', rownames(spatial), value = TRUE)
  }

  results.collector[["mito.genes"]] <- mito.genes
  results.collector[["ribo.genes"]] <- ribo.genes
  results.collector[["ig.genes"]] <- ig.genes

  mito.percent <- Matrix::colSums(spatial[mito.genes, ]) / Matrix::colSums(spatial)
  mito.percent[is.na(mito.percent)] = 0
  ribo.percent <- Matrix::colSums(spatial[ribo.genes, ]) / Matrix::colSums(spatial)
  ribo.percent[is.na(ribo.percent)] = 0
  ig.percent <- Matrix::colSums(spatial[ig.genes, ]) / Matrix::colSums(spatial)
  ig.percent[is.na(ig.percent)] = 0

  results.collector[["mito.percent"]] <- mito.percent
  results.collector[["ribo.percent"]] <- ribo.percent
  results.collector[["ig.percent"]] <- ig.percent

  spatial <- Seurat::AddMetaData(spatial, mito.percent, col.name = "mito.percent")
  spatial <- Seurat::AddMetaData(spatial, ribo.percent, col.name = "ribo.percent")
  spatial <- Seurat::AddMetaData(spatial, ig.percent, col.name = "ig.percent")

  p0 <- Spatial_Plot(spatial,
                     feature = NULL,
                     show.tissue = T,
                     crop = crop,
                     title = sampleName,
                     discrete = T,
                     base.size = 8,
                     ...)

  ## ------ Plot nUMI & nGene ----------
  myPalette <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "Spectral")))

  p.nUMI <- Spatial_Plot(spatial,
                         feature = "nCount_Spatial",
                         show.tissue = F,
                         crop = crop,
                         title = sampleName,
                         legend.title = "nUMI",
                         colors = myPalette(100),
                         discrete = F,
                         base.size = 8,
                         pt.size = 1.8,
                         ...)

  plots.collector[["nUMI"]] <- p.nUMI

  p.nGene <- Spatial_Plot(spatial,
                          feature = "nFeature_Spatial",
                          show.tissue = F,
                          crop = crop,
                          title = sampleName,
                          legend.title = "nGene",
                          colors = myPalette(100),
                          discrete = F,
                          base.size = 8,
                          pt.size = 1.8,
                          ...)

  plots.collector[["nGene"]] <- p.nGene

  p.c <- cbind(ggplot2::ggplotGrob(p0),
               ggplot2::ggplotGrob(p.nUMI),
               ggplot2::ggplotGrob(p.nGene),
               size = "last")

  plots.collector[["basic"]] <- p.c

  ggplot2::ggsave(filename = file.path(savePath_qc, "basic.png"),
                  p.c, height = 5,
                  width = 16,
                  dpi = 300)

  p.nUMI.hist <- histPlot(spatial, feature = "nCount_Spatial", color = "#d56f6d", xlines = c()) +
    labs(title = sampleName)
  p.nGene.hist <- histPlot(spatial, feature = "nFeature_Spatial", color = "#65a55d", xlines = c()) +
    labs(title = sampleName)
  p.c <- cbind(ggplotGrob(p.nUMI.hist), ggplotGrob(p.nGene.hist), size = "last")

  plots.collector[["nUMI_hist"]] <- p.nUMI.hist
  plots.collector[["nGene_hist"]] <- p.nGene.hist
  plots.collector[["basic_hist"]] <- p.c

  ggplot2::ggsave(filename = file.path(savePath_qc, "nUMI-nGene-hist-npc.png"),
                  p.c,
                  dpi = 300,
                  height = 2.5,
                  width = 8)

  p.mito.hist <- histPlot(spatial,
                          feature = "mito.percent",
                          color = "#d56f6d",
                          xlines = c()) +
    labs(title = sampleName)
  p.mito.st <- Spatial_Plot(spatial,
                            feature = "mito.percent",
                            show.tissue = F,
                            crop = crop,
                            title = sampleName,
                            legend.title = "mito.percent",
                            colors = myPalette(100),
                            discrete = F,
                            base.size = 8,
                            pt.size = 1.6,
                            ...)
  p.mito <- plot_grid(p.mito.hist, p.mito.st)

  plots.collector[["mito_hist"]] <- p.mito.hist
  plots.collector[["mito_st"]] <- p.mito.st
  plots.collector[["mito_comb"]] <- p.mito

  ggplot2::ggsave(filename = file.path(savePath_qc, "mito-hist.png"),
                  p.mito.hist,
                  dpi = 300,
                  height = 5,
                  width = 6)
  ggplot2::ggsave(filename = file.path(savePath_qc, "mito-sp.png"),
                  p.mito.st,
                  dpi = 300,
                  height = 5,
                  width = 6)
  ggplot2::ggsave(filename = file.path(savePath_qc, "mito-npc.png"),
                  p.mito,
                  dpi = 300,
                  height = 5,
                  width = 11)

  p.ribo.hist <- histPlot(spatial,
                          feature = "ribo.percent",
                          color = "#65a55d",
                          xlines = c()) +
    labs(title = sampleName)
  p.ribo.st <- Spatial_Plot(spatial,
                            feature = "ribo.percent",
                            show.tissue = F,
                            crop = crop,
                            title = sampleName,
                            legend.title = "ribo.percent",
                            colors = myPalette(100),
                            discrete = F,
                            base.size = 8,
                            pt.size = 1.6,
                            ...)
  p.ribo <- plot_grid(p.ribo.hist, p.ribo.st)

  plots.collector[["ribo_hist"]] <- p.ribo.hist
  plots.collector[["ribo_st"]] <- p.ribo.st
  plots.collector[["ribo_comb"]] <- p.ribo

  ggplot2::ggsave(filename = file.path(savePath_qc, "ribo-hist.png"),
                  p.ribo.hist,
                  dpi = 300,
                  height = 5,
                  width = 6)
  ggplot2::ggsave(filename = file.path(savePath_qc, "ribo-sp.png"),
                  p.ribo.st,
                  dpi = 300,
                  height = 5,
                  width = 6)
  ggplot2::ggsave(filename = file.path(savePath_qc, "ribo-npc.png"),
                  p.ribo,
                  dpi = 300,
                  height = 5,
                  width = 11)

  p.ig.hist <- histPlot(spatial,
                        feature = "ig.percent",
                        color = "#0099CC",
                        xlines = c()) +
    labs(title = sampleName)
  p.ig.st <- Spatial_Plot(spatial,
                          feature = "ig.percent",
                          show.tissue = F,
                          crop = crop,
                          title = sampleName,
                          legend.title = "ig.percent",
                          colors = myPalette(100),
                          discrete = F,
                          base.size = 8,
                          pt.size = 1.6,
                          ...)
  p.ig <- plot_grid(p.ig.hist, p.ig.st)

  plots.collector[["ig_hist"]] <- p.ig.hist
  plots.collector[["ig_st"]] <- p.ig.st
  plots.collector[["ig_comb"]] <- p.ig

  ggplot2::ggsave(filename = file.path(savePath_qc, "ig-hist.png"),
                  p.ig.hist,
                  dpi = 300,
                  height = 5,
                  width = 6)
  ggplot2::ggsave(filename = file.path(savePath_qc, "ig-sp.png"),
                  p.ig.st,
                  dpi = 300,
                  height = 5,
                  width = 6)
  ggplot2::ggsave(filename = file.path(savePath_qc, "ig-npc.png"),
                  p.ig,
                  dpi = 300,
                  height = 5,
                  width = 11)


  ## ------ Gene statistic ------
  files <- list.files(file.path(dataPath, "filtered_feature_bc_matrix"))
  if(!dir.exists(file.path(savePath_data, "filtered_feature_bc_matrix"))){
    dir.create(file.path(savePath_data, "filtered_feature_bc_matrix"), recursive = T)
  }
  for(i in 1 : length(files)){
    file.copy(from = file.path(file.path(dataPath, "filtered_feature_bc_matrix"), files[[i]]),
              to = file.path(file.path(savePath_data, "filtered_feature_bc_matrix"), files[[i]]))
  }

  gene.manifest <- getGeneManifest(file.path(dataPath, "filtered_feature_bc_matrix"))
  gene.manifest$Symbol <- make.unique(names = gene.manifest$Symbol)
  gene.manifest <- addGeneAnno(gene.manifest = gene.manifest,
                               species = species)

  expr.data <- as.data.frame(t(Seurat::Read10X(file.path(dataPath, "filtered_feature_bc_matrix"))))

  nSpot <- Matrix::colSums(expr.data > 0)
  detect.rate <- nSpot / dim(expr.data)[1]

  expr.frac <- expr.data / Matrix::rowSums(expr.data)
  prop.median <- rep(0, length(nSpot))
  names(prop.median) <- names(nSpot)
  tmp.sel.genes <- names(nSpot)[nSpot >= dim(expr.frac)[1] / 2]
  prop.median[tmp.sel.genes] <- apply(expr.frac[, tmp.sel.genes], 2, stats::median)

  gene.manifest$nSpot <- nSpot
  gene.manifest$detect.rate <- detect.rate
  gene.manifest$prop.median <- prop.median

  p.geneProp <- genePropPlot(gene.manifest, t(expr.frac))

  plots.collector[["geneProp"]] <- p.geneProp

  ggplot2::ggsave(filename = file.path(savePath_qc, "geneProp.png"),
                  p.geneProp,
                  dpi = 500,
                  height = 8,
                  width = 8)


  gene.colors <- c(other = "#d56f6d", ribosome = "#f29721", dissociation = "#65a55d", mitochondrial = "#3778bf")
  gene.manifest$Annotation <- factor(gene.manifest$Annotation, ordered = T,
                                     levels = c("mitochondrial", "ribosome", "dissociation", "other"))
  p.gene.rel <- ggplot2::ggplot() +
    ggplot2::geom_point(gene.manifest,
                        mapping = aes(x = detect.rate, y = prop.median, color = Annotation),
                        alpha = 0.3) +
    ggplot2::scale_color_manual(values = gene.colors,
                                name = "Gene") +
    ggplot2::coord_fixed(1/4) +
    ggplot2::scale_y_continuous(trans = 'log10',
                                limits = c(0.00001, 0.1),
                                breaks = c(0.00001, 0.0001, 0.001, 0.01, 0.1),
                                labels = c("0.00001", "0.0001", "0.001", "0.01", "0.1")) +
    ggplot_config(base.size = 6) +
    ggplot2::theme(legend.justification = c(0,1),
                   legend.position = c(0.01,1)) +
    ggplot2::guides(color = guide_legend(override.aes = list(size = 2),
                                         keywidth = 0.1,
                                         keyheight = 0.15,
                                         default.unit = "inch")) +
    ggplot2::labs(x = "Gene detection rate in cells", y = "Median of gene proportion in cells")

  plots.collector[["gene_dr_med"]] <- p.gene.rel

  ggplot2::ggsave(filename =  file.path(savePath_qc, "gene-dr-med.png"),
                  p.gene.rel,
                  dpi = 300,
                  height = 5,
                  width = 5)


  p.nSpot <-
    histPlot.data.frame(gene.manifest,
                        value = "nSpot",
                        color = "#B3B3B3",
                        xlines = c()) +
    ggplot2::scale_y_continuous(trans = 'log10') +
    ggplot2::labs(y = "Gene number")

  plots.collector[["nSpot"]] <- p.nSpot

  ggplot2::ggsave(filename = file.path(savePath_qc, "nSpot-hist.png"),
                  p.nSpot,
                  dpi = 300,
                  height = 5,
                  width = 8)

  results.collector[["gene.manifest"]] <- gene.manifest

  write.table(gene.manifest,
              file = file.path(savePath_data, "geneManifest.txt"),
              quote = F,
              sep = "\t",
              row.names = F)

  saveRDS(spatial, file.path(savePath_data, "spatial_stat.RDS"))

  if(genReport){
    knitr::knit(file.path(system.file(package = "stCancer"), "rmd/stStat.Rmd"),
                file.path(savePath_basic, "report-stStat.md"))
    markdown::markdownToHTML(file.path(savePath_basic, 'report-stStat.md'),
                             file.path(savePath_basic, 'report-stStat.html'))
  }

  return(list(object = spatial,
              results = results.collector,
              plots = plots.collector))
}


#' Load_Spatial
#'
#' Load a 10x Genomics Visium Spatial Experiment into a Seurat object (similar to Load10X_Spatial in Seurat)
#'
#' @param data.dir Directory containing the H5 file specified by filename and the image data in a subdirectory called spatial.
#' Under this path, folders 'filtered_feature_bc_matrix' and 'raw_feature_bc_matrix' exist generally.
#' @param assay Name of the initial assay
#' @param slice Name for the stored image of the tissue slice.
#' @param sampleName A character string giving a label for this sample.
#' @param rm.isolated A logical value indicating whether to remove isolated spots.
#' @param region.threshold A integer value, if the area of connected domain of spots are less than `region.threshold`,
#' they will be removed.
#' @param filter.matrix Only keep spots that have been determined to be over tissue.
#' @param to.upper Converts all feature names to upper case. Can be useful when analyses require comparisons between human and mouse gene names for example.
#' @param h5 A logical value indicating whether to read h5 files, whose contents are the same as 'filtered_feature_bc_matrix' or 'raw_feature_bc_matrix'
#' @param ...
#'
#'
#' @return A Seurat object
#' @export
#'
#' @import Seurat
#'
Load_Spatial <- function(data.dir,
                         assay = "Spatial",
                         slice = "slice1",
                         sampleName = "SeuratProject",
                         rm.isolated = TRUE,
                         region.threshold = 3,
                         h5 = FALSE,
                         filter.matrix = TRUE,
                         to.upper = FALSE,
                         ...){
  if(length(data.dir) > 1){
    warning("'Load_Spatial' accepts only one 'data.dir'", immediate. = TRUE)
    data.dir <- data.dir[1]
  }
  if(!h5){
    if(filter.matrix){
      data <- as.data.frame(Seurat::Read10X(file.path(data.dir, "filtered_feature_bc_matrix/")))
    }else{
      data <- as.data.frame(Seurat::Read10X(file.path(data.dir, "raw_feature_bc_matrix/")))
    }
  }else{
    if(filter.matrix){
      data <- as.data.frame(Seurat::Read10X_h5(file.path(data.dir, "filtered_feature_bc_matrix.h5")))
    }else{
      data <- as.data.frame(Seurat::Read10X_h5(file.path(data.dir, "raw_feature_bc_matrix.h5")))
    }
  }
  if(to.upper){
    rownames(data) <- toupper(rownames(data))
  }
  t.list <- Read_Image(image.dir = file.path(data.dir, "spatial"),
                       filter.matrix = filter.matrix,
                       rm.isolated = rm.isolated,
                       region.threshold = region.threshold)
  image <- t.list[[1]]

  tissue.positions <- read.csv(
    file = file.path(file.path(data.dir, "spatial"), 'tissue_positions_list.csv'),
    col.names = c('barcodes', 'tissue', 'row', 'col', 'imagerow', 'imagecol'),
    header = FALSE,
    as.is = TRUE,
    row.names = 1
  )
  image@scale.factors$crop.row.min <- min(tissue.positions$imagerow)
  image@scale.factors$crop.row.max <- max(tissue.positions$imagerow)
  image@scale.factors$crop.col.min <- min(tissue.positions$imagecol)
  image@scale.factors$crop.col.max <- max(tissue.positions$imagecol)

  data <- data[, t.list[[2]]]
  object <- Seurat::CreateSeuratObject(counts = data,
                                       assay = assay)
  image <- image[Seurat::Cells(object)]
  Seurat::DefaultAssay(image) <- assay
  object[[slice]] <- image
  object@project.name <- sampleName
  return(object)
}


#' Read_Image
#'
#' Load a 10x Genomics Visium Image (similar to Read10X_Image in Seurat)
#'
#' @param image.dir Path to directory with 10X Genomics visium image data; should include files tissue_lowres_iamge.png, scalefactors_json.json and tissue_positions_list.csv
#' @param filter.matrix Filter spot/feature matrix to only include spots that have been determined to be over tissue.
#' @param rm.isolated A logical value indicating whether to remove isolated spots.
#' @param region.threshold A integer value, if the area of connected domain of spots are less than `region.threshold`,
#' they will be removed.
#' @param ...
#'
#'
#' @return A VisiumV1 object
#' @export
#'
#' @import Seurat
#' @importFrom png readPNG
#' @importFrom jsonlite fromJSON
#'
Read_Image <- function(image.dir,
                       filter.matrix = TRUE,
                       rm.isolated = T,
                       region.threshold = 3,
                       ...){
  image <- png::readPNG(source = file.path(image.dir, 'tissue_lowres_image.png'))
  scale.factors <- jsonlite::fromJSON(txt = file.path(image.dir, 'scalefactors_json.json'))
  tissue.positions <- read.csv(
    file = file.path(image.dir, 'tissue_positions_list.csv'),
    col.names = c('barcodes', 'tissue', 'row', 'col', 'imagerow', 'imagecol'),
    header = FALSE,
    as.is = TRUE,
    row.names = 1
  )
  if(rm.isolated){
    tissue.positions <- rm_isolated_spot(tissue.positions,
                                         region.threshold = region.threshold)
  }
  if(filter.matrix){
    tissue.positions <- tissue.positions[which(x = tissue.positions$tissue == 1), , drop = FALSE]
  }

  unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
  spot.radius <-  unnormalized.radius / max(dim(x = image))
  return(list(new(
    Class = 'VisiumV1',
    image = image,
    scale.factors = scalefactors(
      spot = scale.factors$tissue_hires_scalef,
      fiducial = scale.factors$fiducial_diameter_fullres,
      hires = scale.factors$tissue_hires_scalef,
      lowres = scale.factors$tissue_lowres_scalef
    ),
    coordinates = tissue.positions,
    spot.radius = spot.radius
  ),
  rownames(tissue.positions)))
}


#' rm_isolated_spot
#'
#' Remove isolated spots
#'
#' @param spm the data.frame of tissue.possition
#' @param region.threshold A integer value, if the area of connected domain of spots are less than `region.threshold`,
#' they will be removed.
#'
#' @return the data.frame of tissue.possition after removing isolated areas
#'
#' @import Matrix
#'
rm_isolated_spot <- function(spm,
                             region.threshold = 3){
  t_spm <- spm[spm$tissue == 1, ]
  # Connected domain labeling
  t_spm$label <- -1
  label_id <- 1
  for(j in 1 : nrow(t_spm)){
    if(t_spm[j, "label"] != -1){
      next
    }

    t_spm[j, "label"] <- label_id

    j_row <- t_spm[j, ]$row
    j_col <- t_spm[j, ]$col
    j_barcode <- rownames(t_spm[j, ])
    nei <- get_neighbors_indices(j_row, j_col)
    # area list
    spot_list <- list()
    spot_list[[1]] <- c(j_row, j_col)
    while(length(spot_list) != 0){
      ii_spot <- spot_list[[1]]
      spot_list[[1]] <- NULL
      nei <- get_neighbors_indices(ii_spot[1], ii_spot[2])
      for(k in 1 : length(nei)){
        i_spm <- t_spm[t_spm$row == nei[[k]][1] &
                         t_spm$col == nei[[k]][2], ]
        if(nrow(na.omit(i_spm)) == 0){
          next
        }
        if(i_spm$label == -1){
          spot_list[[length(spot_list) + 1]] <- c(i_spm$row, i_spm$col)
          t_spm[t_spm$row == nei[[k]][1] &
                  t_spm$col == nei[[k]][2], "label"] <- label_id
        }
      }
    }
    label_id <- label_id + 1
  }
  # judge the size of area
  for(j in 1 : label_id){
    i_spm <- t_spm %>% subset(label == j)
    i_spm <- na.omit(i_spm)
    if(nrow(i_spm) <= region.threshold & nrow(i_spm) != 0){
      i_barcode <- rownames(i_spm)
      spm[i_barcode, ]$tissue <- 0
    }
  }

  return(spm)
}
