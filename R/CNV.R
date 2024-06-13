#' runCNVAnalysis
#'
#' CNV analysis
#'
#' @param SeuratObject A Seurat object
#' @param savePath A path to save the results files
#' @param analysis.func The function used to analyze CNV
#' @param ref.object A seuratObject used as the reference if using InferCNV
#' @param species "human" or "mouse", only human data can be used in CNV analysis by copyKAT
#'
#'
#' @return A list containing a seuratObject, a list of results of intermediate process,
#' a list of plots and a logical value to show whether the process completed successfully
#' @export
#'
#' @import Seurat ggplot2 ComplexHeatmap png
#' @importFrom dplyr "%>%"
#' @importFrom circlize colorRamp2
runCNVAnalysis <- function(SeuratObject,
                           savePath,
                           analysis.func = "copyKAT",
                           ref.object = NULL,
                           species = "human",
                           ...){

    # collect results and plots
    results.collector <- list()
    plots.collector <- list()
    bool.completed <- FALSE

    if(!((species == "human") | (species == "mouse" & analysis.func == "inferCNV"))){
        warning(paste0("The species of ", species, " is not supported for interaction analysis."))
        return(list(object = SeuratObject,
                    results = results.collector,
                    plots = plots.collector,
                    bool.completed = bool.completed))
    }

    if(analysis.func == "inferCNV"){
        return(runInferCNV(SeuratObject,
                           savePath,
                           ref.object = ref.object,
                           species = species,
                           ...))
    }else if(analysis.func == "copyKAT"){
        return(results <- runCopyKAT(SeuratObject,
                                     savePath,
                                     ...))
    }else{
        stop(paste0("The function ", analysis.func, " is not supported. Please choose copyKAT or inferCNV."))
        return(1)
    }
}


#' runInferCNV
#'
#' Analyze CNV using the method of InferCNV
#'
#' @param SeuratObject A seuratObejct to analyze CNV
#' @param savePath A path to store the results
#' @param cluster.info The label to store cluster information in `SeuratObject@meta.data`
#' @param ref.object A seuratObject used as the reference if using InferCNV
#' @param species "human" or "mouse"
#' @param genome
#'
#' @return A list containing a seuratObject, a list of results of intermediate process,
#' a list of plots and a logical value to show whether the process completed successfully
#' @export
#'
runInferCNV <- function(SeuratObject,
                        savePath,
                        ref.object,
                        cluster.info = "seurat_clusters",
                        species = "human",
                        genome = "hg19"){

    # collect results and plots
    results.collector <- list()
    plots.collector <- list()
    bool.completed <- FALSE

    sampleName <- SeuratObject@project.name

    if(is.null(ref.object)){
        stop("The reference is necessary!")
        return(1)
    }

    ori.cluster <- SeuratObject@meta.data[[cluster.info]]
    clusters <- unique(ori.cluster)
    clusters <- as.numeric(as.character(clusters))

    ref.object <- SCTransform(ref.object, assay = ref.object@active.assay, verbose = FALSE)
    ref.object <- RunPCA(ref.object, assay = ref.object@active.assay, verbose = FALSE)
    ref.object <- FindVariableFeatures(ref.object, verbose = FALSE)
    ref.object <- FindNeighbors(ref.object, reduction = "pca", dims = 1:30, verbose = FALSE)
    ref.object <- FindClusters(ref.object, resolution = 0.5, verbose = FALSE)

    t.results <- runMalignancy(SeuratObject,
                               savePath = savePath,
                               assay = "Spatial",
                               cutoff = 0.1,
                               minspot = 3,
                               p.value.cutoff = 0.05,
                               coor.names = c("tSNE_1", "tSNE_2"),
                               ref.object = ref.object,
                               species = species,
                               genome = genome)

    results.collector[["inferCNV"]] <- t.results

    # Draw
    gene.chr <- read.table(system.file("txt", "gene-chr-hg38.txt", package = "stCancer"),
                           col.names = c("EnsemblID", "CHR", "C_START", "C_STOP"),
                           stringsAsFactors = F)

    chr.colors <- c("#cb7c77", "#68d359", "#6a7dc9", "#c9d73d", "#c555cb",
                    "#d7652d", "#7cd5c8", "#c49a3f", "#f6a97a", "#507d41",
                    "#c8b693", "#90353b", "#aed688", "#674c2a", "#c6a5cc",
                    "#1B9E77", "#c5383c", "#0081d1", "#d1b900", "#502e71",
                    "#cf4c8b", "#798234", "#6b42c8", "#5d8d9c", "#666666")

    cnv.data <- read.table(file.path(savePath, "inferCNV-observation.txt"),
                           header = T,
                           stringsAsFactors = F)
    colnames(cnv.data) <- sapply(colnames(cnv.data), function(x){
        return(gsub("\\.", "-", x))
    })

    spot.anno <- data.frame(Cluster = ori.cluster,
                            row.names = colnames(malig.object))

    gene.anno <- subset(gene.chr, EnsemblID %in% rownames(cnv.data))
    cnv.data <- cnv.data[gene.anno$EnsemblID, ]
    rownames(gene.anno) <- gene.anno$EnsemblID
    gene.anno <- gene.anno[, c("CHR"), drop = F]
    gene.anno$CHR <- factor(gene.anno$CHR, levels = unique(gene.anno$CHR))

    cur.color <- malig.object@meta.data[["cluster.color"]]
    names(cur.color) <- ori.cluster

    chr.color <- chr.colors[1:length(unique(gene.anno$CHR))]
    names(chr.color) <- unique(gene.anno$CHR)

    p.cnvdata <- cnv.data[, colnames(malig.object)]
    up.lim <- quantile(as.matrix(p.cnvdata), 0.95)
    p.cnvdata[p.cnvdata > up.lim] <- up.lim
    p.colors <- circlize::colorRamp2(c(0.9, 0.98, 1.02, 1.1), c("blue", "white", "white", "red"))
    spot.color <- list(Cluster = cur.color)
    gene.color <- list(CHR = chr.color)
    spot.ha <- ComplexHeatmap::rowAnnotation(df = spot.anno,
                                             col = spot.color,
                                             show_annotation_name = F)
    gene.ha <- ComplexHeatmap::HeatmapAnnotation(df = gene.anno,
                                                 col = gene.color,
                                                 annotation_name_side = "left")

    png(filename = file.path(savePath, "cnv-heat.png"),
        width = 3000, height = 1500, res = 300)
    suppressWarnings(
        p <- ComplexHeatmap::Heatmap(t(p.cnvdata),
                                     column_title = sampleName,
                                     show_column_names = F,
                                     show_row_names = F,
                                     cluster_columns = F,
                                     show_row_dend = F,
                                     name = "CNV",
                                     col = p.colors,
                                     left_annotation = spot.ha,
                                     top_annotation = gene.ha,
                                     row_split = spot.anno$Cluster,
                                     use_raster = TRUE)
    )
    print(p)
    garbage <- dev.off()

    plots.collector[["CNV_heat"]] <- p

    bool.completed <- TRUE
    return(list(object = object,
                results = results.collector,
                plots = plots.collector,
                bool.completed = bool.completed))
}


#' runCopyKAT
#'
#' Analyze CNV using R package copyKAT
#' Gao, R., Bai, S., Henderson, Y. C., Lin, Y., Schalck, A., Yan, Y., Kumar,
#' T., Hu, M., Sei, E., Davis, A., Wang, F., Shaitelman, S. F., Wang, J. R.,
#' Chen, K., Moulder, S., Lai, S. Y. & Navin, N. E. (2021).
#' Delineating copy number and clonal substructure in human tumors from
#' single-cell transcriptomes. Nat Biotechnol. doi:10.1038/s41587-020-00795-2.
#'
#' @param SeuratObject A seuratObejct to analyze CNV
#' @param savePath A path to store the results
#' @param cluster.info The label to store cluster information in `SeuratObject@meta.data`
#' @param id.type The parameter of copyKAT
#' @param ngene.chr The parameter of copyKAT
#' @param win.size The parameter of copyKAT
#' @param KS.cut The parameter of copyKAT
#' @param n.cores The parameter of copyKAT
#' @param sam.name The parameter of copyKAT
#' @param distance The parameter of copyKAT
#' @param norm.cell.names The parameter of copyKAT
#' @param output.set The parameter of copyKAT
#'
#' @return A list containing a seuratObject, a list of results of intermediate process,
#' a list of plots and a logical value to show whether the process completed successfully
#' @export
#'
#' @import copykat
runCopyKAT <- function(SeuratObject,
                       savePath,
                       cluster.info = "seurat_clusters",
                       id.type="S",
                       ngene.chr=5,
                       win.size=25,
                       KS.cut=0.1,
                       n.cores = 1,
                       sam.name="test",
                       distance="euclidean",
                       norm.cell.names="",
                       output.seg="FLASE",
                       ...){

    # collect results and plots
    results.collector <- list()
    plots.collector <- list()
    bool.completed <- FALSE

    sampleName <- SeuratObject@project.name

    full.anno <<- copykat::full.anno
    cyclegenes <<- copykat::cyclegenes
    DNA.hg20 <<- copykat::DNA.hg20
    exp.rawdata <<- copykat::exp.rawdata

    tmp.path <- getwd()
    setwd(file.path(savePath))

    if(class(SeuratObject@assays[[1]])[1] == 'Assay5'){
        rawmat <- Seurat::GetAssayData(SeuratObject, assay = SeuratObject@active.assay, layer = 'counts')
    }else{
        rawmat <- Seurat::GetAssayData(SeuratObject, assay = SeuratObject@active.assay, slot = 'counts')
    }

    copykat.test <- copykat::copykat(rawmat=rawmat,
                                     id.type=id.type,
                                     ngene.chr=ngene.chr,
                                     win.size=win.size,
                                     KS.cut=KS.cut,
                                     n.cores = n.cores,
                                     sam.name=sam.name,
                                     distance=distance,
                                     norm.cell.names=norm.cell.names,
                                     output.seg=output.seg,
                                     ...)

    saveRDS(copykat.test, file.path(savePath, 'copykat.rds'))
    # rm(full.anno, cyclegenes, DNA.hg20, exp.rawdata)
    setwd(tmp.path)

    results.collector[["copyKAT"]] <- copykat.test

    pred.test <- data.frame(copykat.test$prediction)
    CNA.test <- data.frame(copykat.test$CNAmat)

    # Draw
    gene.chr <- read.table(system.file("txt", "gene-chr-hg38.txt", package = "stCancer"),
                           col.names = c("EnsemblID", "CHR", "C_START", "C_STOP"),
                           stringsAsFactors = F)

    chr.colors <- c("#cb7c77", "#68d359", "#6a7dc9", "#c9d73d", "#c555cb",
                    "#d7652d", "#7cd5c8", "#c49a3f", "#f6a97a", "#507d41",
                    "#c8b693", "#90353b", "#aed688", "#674c2a", "#c6a5cc",
                    "#1B9E77", "#c5383c", "#0081d1", "#d1b900", "#502e71",
                    "#cf4c8b", "#798234", "#6b42c8", "#5d8d9c", "#666666")

    pred.tag <- data.frame(malig.pred = rep("unknown", ncol(SeuratObject)),
                           row.names = colnames(SeuratObject))
    pred.tag[rownames(pred.test), ] <- pred.test$copykat.pred

    p.cnvdata <- CNA.test[, -c(1, 2, 3)]
    p.cnvdata <- p.cnvdata + 1
    rownames(p.cnvdata) <- CNA.test[, 3]
    colnames(p.cnvdata) <- sapply(colnames(p.cnvdata), function(x){
        return(gsub("\\.", "-", x))
    })

    subObject <- SeuratObject[, colnames(p.cnvdata)]
    cur.color <- subObject@meta.data[["cluster.color"]]
    names(cur.color) <- subObject@meta.data[[cluster.info]]
    spot.color <- list(Cluster = cur.color)
    spot.anno <- data.frame(Cluster = subObject@meta.data[[cluster.info]],
                            row.names = colnames(p.cnvdata))

    gene.anno <- data.frame(CHR = CNA.test[, 1],
                            row.names = CNA.test[, 3])
    gene.anno$CHR <- sapply(gene.anno$CHR, function(x){
        return(paste0("chr", x))
    })

    chr.color <- chr.colors[1:length(unique(gene.anno$CHR))]
    names(chr.color) <- unique(gene.anno$CHR)
    gene.color <- list(CHR = chr.color)

    up.lim <- quantile(as.matrix(p.cnvdata), 0.95)
    p.cnvdata[p.cnvdata > up.lim] <- up.lim
    p.colors <- circlize::colorRamp2(c(0.9, 0.98, 1.02, 1.1), c("blue", "white", "white", "red"))

    spot.ha <- ComplexHeatmap::rowAnnotation(df = spot.anno,
                                             col = spot.color,
                                             show_annotation_name = F)
    gene.ha <- ComplexHeatmap::HeatmapAnnotation(df = gene.anno,
                                                 col = gene.color,
                                                 annotation_name_side = "left")

    png(filename = file.path(savePath, "cnv-heat.png"),
        width = 3000, height = 1500, res = 300)
    suppressWarnings(
        p <- ComplexHeatmap::Heatmap(t(p.cnvdata),
                                     column_title = sampleName,
                                     show_column_names = F,
                                     show_row_names = F,
                                     cluster_columns = F,
                                     show_row_dend = F,
                                     name = "CNV",
                                     col = p.colors,
                                     left_annotation = spot.ha,
                                     top_annotation = gene.ha,
                                     row_split = spot.anno$Cluster,
                                     use_raster = TRUE)
    )
    print(p)
    garbage <- dev.off()

    plots.collector[["CNV_heat"]] <- p

    # label
    SeuratObject <- AddMetaData(SeuratObject, pred.tag, col.name = "malig.pred")

    colors <- c("#808080", # gray
                "#9370DB", # purple
                "#B0E0E6") # blue
    names(colors) <- c("unknown", "aneuploid", "diploid")

    malig.plot <- SpatialDim_Plot(SeuratObject,
                                  group.by = 'malig.pred',
                                  crop = T,
                                  legend.color = colors)

    # malig.plot <- Spatial_Plot(SeuratObject,
    #                            "malig.pred",
    #                            discrete = T,
    #                            show.tissue = T,
    #                            crop = T,
    #                            pt.size = 2.0,
    #                            colors = colors)

    plots.collector[["malig.pred"]] <- malig.plot

    ggsave(file.path(savePath, "cnv-sp.png"),
           malig.plot,
           height = 5,
           width = 7)

    bool.completed <- TRUE
    return(list(object = object,
                results = results.collector,
                plots = plots.collector,
                bool.completed = bool.completed))
}






## functions related to inferCNV
#' prepareCNV
#'
#' Get common genes in 'expr.data' and 'ref.data'
#'
#' @param SeuratObject A seurat object
#' @param gene.manifest A dataframe of gene information
#' @param assay Assay
#' @param ref.data The reference seurat object data
#' @param species Species
#' @param genome Genome, hg38, hg19 or mm10
#'
#'
#' @return A list
#' @export
#'
#' @import Seurat
#' @importFrom dplyr distinct
prepareCNV <- function(SeuratObject,
                       gene.manifest,
                       assay = "Spatial",
                       ref.object = NULL,
                       species = "human",
                       genome = "hg19"){
    ## gene.chr
    if(species == "human"){
        if(genome == "hg38"){
            gene.chr <- read.table(system.file("txt", "gene-chr-hg38.txt", package = "stCancer"),
                                   col.names = c("EnsemblID", "CHR", "C_START", "C_STOP"),
                                   stringsAsFactors = F)
        }else if(genome == "hg19"){
            gene.chr <- read.table(system.file("txt", "gene-chr-hg19.txt", package = "stCancer"),
                                   col.names = c("EnsemblID", "CHR", "C_START", "C_STOP"),
                                   stringsAsFactors = F)
        }else{
            stop("Error in 'runInferCNV': ", genome, " is not allowed for 'genome'.\n")
        }
    }else if(species == "mouse"){
        if(genome == "mm10"){
            gene.chr <- read.table(system.file("txt", "gene-chr-mm10.txt", package = "stCancer"),
                                   col.names = c("EnsemblID", "CHR", "C_START", "C_STOP"),
                                   stringsAsFactors = F)
        }else{
            stop("Error in 'runInferCNV': ", genome, " is not allowed for 'genome'.\n")
        }
    }else{
        stop("Error in 'runInferCNV': ", species, " is not allowed for 'species'.\n")
    }

    ## reference.data
    if(is.null(ref.object)){
        if(species == "human"){
            ref.data <- readRDS(system.file("rds", "cnvRef_Data-HM.RDS", package = "stCancer"))
        }else if(species == "mouse"){
            ref.data <- readRDS(system.file("rds", "cnvRef_Data-boneMarrow-MS.RDS", package = "stCancer"))
        }
    }else{
        ref.data <- GetAssayData(ref.object, assay = assay, slot = "counts")
    }

    expr.data <- GetAssayData(SeuratObject, assay = assay, slot = "counts")

    ref.anno <- data.frame(barcode = colnames(ref.data),
                           spotAnno = "Reference",
                           stringsAsFactors = F)
    rownames(ref.anno) <- ref.anno$barcode

    ## combine data
    com.genes <- intersect(rownames(expr.data), rownames(ref.data))
    ref.data <- ref.data[com.genes, ]
    expr.data <- expr.data[com.genes, ]

    gene.manifest <- gene.manifest %>% distinct(Symbol, .keep_all = T)

    rownames(gene.manifest) <- gene.manifest$Symbol
    rownames(expr.data) <- gene.manifest[rownames(expr.data), ]$EnsemblID

    ## combine spot.anno
    spot.anno <- data.frame(barcode = colnames(expr.data),
                            spotAnno = "Experiment",
                            stringsAsFactors = F)
    rownames(spot.anno) <- spot.anno$barcode

    spot.anno <- rbind(spot.anno, ref.anno)

    ## common genes between expr.data and gene.chr
    com.genes <- intersect(rownames(expr.data), gene.chr$EnsemblID)
    gene.chr <- subset(gene.chr, EnsemblID %in% com.genes)

    # gene.chr <- gene.chr[with(gene.chr, order(CHR, C_START, C_STOP)), ]
    expr.data <- cbind(as.matrix(expr.data), ref.data)

    expr.data <- expr.data[gene.chr$EnsemblID, ]
    rownames(gene.chr) <- gene.chr$EnsemblID

    return(list(expr.data = expr.data,
                gene.chr = gene.chr,
                spot.anno = spot.anno))
}


# Remove genes whose mean expression are less than 'cutoff' and total counts are less than 'minspot'
rmGeneForCNV <- function(cnvList, cutoff = 0.1, minspot = 3){
    gene.mean <- Matrix::rowMeans(cnvList$expr.data)
    gene.sum <- Matrix::rowSums(cnvList$expr.data > 0)
    genes.sel <- rownames(cnvList$expr.data)[gene.mean >= cutoff & gene.sum >= minspot]

    cnvList$expr.data <- cnvList$expr.data[genes.sel, ]
    cnvList$gene.chr <- cnvList$gene.chr[genes.sel, ]

    return(cnvList)
}


normalizeDataForCNV <- function(cnvList){
    expr.data <- cnvList$expr.data

    cs = Matrix::colSums(expr.data)
    expr.data <- t(t(expr.data) / cs)
    normalize_factor <- 10^round(log10(mean(cs)))
    expr.data <- expr.data * normalize_factor

    cnvList$expr.data <- expr.data
    return(cnvList)
}


# Anscombe transform
anscombeTransform <- function(cnvList){
    cnvList$expr.data <- 2 * sqrt(cnvList$expr.data + 3/8)
    return(cnvList)
}



logForCNV <- function(cnvList){
    cnvList$expr.data <- log2(cnvList$expr.data + 1)
    return(cnvList)
}



getAverageBounds <- function(cnvList){
    lower.bound <- mean(apply(cnvList$expr.data, 2, min))
    upper.bound <- mean(apply(cnvList$expr.data, 2, max))
    threshold = mean(abs(c(lower.bound, upper.bound)))
    return(threshold)
}



boundForCNV <- function(cnvList, threshold){
    cnvList$expr.data[cnvList$expr.data > threshold] <- threshold
    cnvList$expr.data[cnvList$expr.data < (-1 * threshold)] <- -1 * threshold
    return(cnvList)
}



smoothOne <- function(ori.data, window.len = window.len){
    half.window <- (window.len - 1) / 2

    pad.data <- c(rep(0, half.window), ori.data, rep(0, half.window))
    bool.data <- c(rep(0, half.window), rep(1, length(ori.data)), rep(0, half.window))

    kernel.vec <- c(1:half.window, half.window + 1, half.window:1)

    sum.data <- stats::filter(pad.data, kernel.vec, sides = 2)
    num.data <- stats::filter(bool.data, kernel.vec, sides = 2)
    sum.data <- sum.data[!is.na(sum.data)]
    num.data <- num.data[!is.na(num.data)]

    smo.data <- sum.data / num.data
    return(smo.data)
}


smoothByChr <- function(cnvList, window.len = 101){
    chrList <- cnvList$gene.chr$CHR
    chrs <- as.character(unique(cnvList$gene.chr$CHR))

    if(window.len < 2){
        cat("- Warning in 'smoothBychr': Window length < 2, returning original data.\n")
        return(cnvList)
    }

    expr.data <- cnvList$expr.data
    for(chr in chrs) {
        # print(chr)
        cur.genes.ix <- which(chrList == chr)
        cur.data <- expr.data[cur.genes.ix, , drop=F]

        if(length(cur.genes.ix) > 1) {
            if(window.len %% 2 == 0){
                window.len = window.len + 1
                cat("- Warning in 'smoothBychr': Window length is even, adding one to 'window.len'.\n")
            }

            smooth.data <- apply(cur.data, 2, smoothOne, window.len = window.len)
            rownames(smooth.data) <- rownames(cur.data)
            colnames(smooth.data) <- colnames(cur.data)

            expr.data[cur.genes.ix, ] <- smooth.data
        }
    }
    cnvList$expr.data <- expr.data
    return(cnvList)
}



centerAcrossChr <- function(cnvList, method = "median"){
    expr.data <- cnvList$expr.data
    if (method == "median") {
        row_median <- apply(expr.data, 2, function(x) { median(x, na.rm=T) } )
        expr.data <- t(apply(expr.data, 1, "-", row_median))
    }
    else {
        row_means <- apply(expr.data, 2, function(x) { mean(x, na.rm=T) } )
        expr.data <- t(apply(expr.data, 1, "-", row_means))
    }
    cnvList$expr.data <- expr.data
    return(cnvList)
}



subtractRefExpr <- function(cnvList, inv_log = TRUE){
    ref.spotNames <- subset(cnvList$spot.anno, spotAnno == "Reference")$barcode

    if (inv_log) {
        ref.means = log2(Matrix::rowMeans(2^cnvList$expr.data[, ref.spotNames] - 1) + 1)
    } else {
        ref.means = Matrix::rowMeans(cnvList$expr.data[, ref.spotNames])
    }

    cnvList$expr.data <- cnvList$expr.data - ref.means
    return(cnvList)
}


invertLog2 <- function(cnvList){
    cnvList$expr.data <- 2^cnvList$expr.data
    return(cnvList)
}



denoiseByRefMeanSd <- function(cnvList, sd_amplifier=1.5){
    expr.data <- cnvList$expr.data
    ref.spotNames <- subset(cnvList$spot.anno, spotAnno == "Reference")$barcode

    mean.ref.vals <- mean(expr.data[, ref.spotNames])
    mean.ref.sd <- mean(apply(expr.data[, ref.spotNames], 2, function(x) sd(x, na.rm=T))) * sd_amplifier

    up.bound <- mean.ref.vals + mean.ref.sd
    low.bound <- mean.ref.vals - mean.ref.sd

    expr.data[expr.data > low.bound & expr.data < up.bound] <- mean.ref.vals

    cnvList$expr.data <- expr.data

    return(cnvList)
}



removeOutliers <- function(cnvList){
    expr.data <- cnvList$expr.data
    up.bound <- mean(apply(expr.data, 2, max))
    low.bound <- mean(apply(expr.data, 2, min))

    expr.data[expr.data < low.bound] <- low.bound
    expr.data[expr.data > up.bound] <- up.bound

    cnvList$expr.data <- expr.data

    return(cnvList)
}



runCNV <- function(SeuratObject,
                   gene.manifest,
                   assay = "Spatial",
                   cutoff = 0.1,
                   minspot = 3,
                   ref.object = NULL,
                   species = "human",
                   genome = "hg19"){

    cnvList <- prepareCNV(SeuratObject = SeuratObject,
                          assay = assay,
                          gene.manifest = gene.manifest,
                          ref.object = ref.object,
                          species = species,
                          genome = genome)
    cnvList <- rmGeneForCNV(cnvList, cutoff = cutoff, minspot = minspot)
    cnvList <- normalizeDataForCNV(cnvList)
    cnvList <- anscombeTransform(cnvList)
    cnvList <- logForCNV(cnvList)


    threshold <- getAverageBounds(cnvList)
    cnvList <- boundForCNV(cnvList, threshold)
    cnvList <- smoothByChr(cnvList, window.len = 101)
    cnvList <- centerAcrossChr(cnvList, method = "median")
    cnvList <- subtractRefExpr(cnvList)
    cnvList <- invertLog2(cnvList)
    cnvList <- denoiseByRefMeanSd(cnvList, sd_amplifier = 1.0)
    cnvList <- removeOutliers(cnvList)

    return(cnvList)
}


getMalignScore <- function(cnvList, spot.type = "Observation", method = "smooth", adjMat = NULL){
    if(spot.type == "Observation"){
        spot.names <- subset(cnvList$spot.anno, spotAnno != "Reference")$barcode
    }else if(spot.type == "Reference"){
        spot.names <- subset(cnvList$spot.anno, spotAnno == "Reference")$barcode
    }

    cur.data <- cnvList$expr.data[, spot.names]

    if(is.null(adjMat) & method == "smooth"){
        cat("- Warning in 'getMalignScore': Adjacent matrix is not provided, and use 'direct' method instead.\n")
        method <- "direct"
    }
    if(method == "smooth"){
        thres <- quantile(adjMat@x, 1- (dim(adjMat)[1] * 10 / length(adjMat@x)))

        indexes <- as.matrix((adjMat > thres) + 0)
        tt <- 0.5 / (rowSums(indexes) - 1)
        tt[is.infinite(tt)] <- 0

        indexes <- indexes * tt
        indexes <- indexes * (1 - diag(rep(1, dim(indexes)[1])))
        diagValue <- rep(0.5, dim(indexes)[1])
        diagValue[tt == 0] <- 1

        indexes <- t(indexes + diag(diagValue))

        new.cur.data <- as.matrix(cur.data) %*% indexes
        malignScore <- colSums((new.cur.data - 1)^2)
        malignScore <- malignScore / dim(new.cur.data)[1]

    }else if(method == "direct"){
        malignScore <- colSums((cur.data - 1)^2)
        malignScore <- malignScore / dim(cur.data)[1]
    }

    names(malignScore) <- colnames(cur.data)

    return(malignScore)
}



malignPlot <- function(obserScore, referScore, malign.thres = NULL){
    scoreDF <- data.frame(malignScore = c(obserScore, referScore),
                          sets = c(rep("Observation", length(obserScore)),
                                   rep("Reference", length(referScore))))
    p <- ggplot() +
        geom_histogram(data = subset(scoreDF, sets == "Observation"),
                       mapping = aes(x = malignScore, fill = "Observation"),
                       bins = 150, alpha = 0.6) +
        geom_histogram(data = subset(scoreDF, sets == "Reference"),
                       mapping = aes(x = malignScore, fill = "Reference"),
                       bins = 150, alpha = 0.6) +
        labs(x = "Malignancy score", y = "Droplets count") +
        scale_fill_manual(name = "spots sets", guide = "legend",
                          values = c("Observation"="#2e68b7", "Reference"="grey")) +
        theme_classic() +
        ggplot_config(base.size = 7) +
        theme(legend.justification = c(1.12,1.12), legend.position = c(1,1))
    if(!is.null(malign.thres)){
        p <- p + geom_vline(xintercept = malign.thres, colour = "red", linetype = "dashed")
    }
    return(p)
}



getBimodalThres <- function(scores){
    x.density <- density(scores)
    d.x.density <- diff(x.density$y)
    d.sign <- (d.x.density > 0) + 0

    ext.pos <- which(d.sign[2:length(d.sign)] - d.sign[1:(length(d.sign)-1)] != 0)
    ext.density <- x.density$y[ext.pos]
    y.max <- max(ext.density)
    if(length(ext.pos) >= 3){
        del.ix <- c()
        for(ei in 2:length(ext.density)){
            if(abs(ext.density[ei] - ext.density[ei - 1]) < y.max * 0.001){
                del.ix <- c(del.ix, ei - 1, ei)
            }
        }
        sel.ix <- !(1:length(ext.density) %in% unique(del.ix))
        ext.density <- ext.density[sel.ix]
        ext.pos <- ext.pos[sel.ix]
    }

    if(length(ext.pos) >= 3){
        t.ext.density <- c(0, ext.density, 0)
        ext.height <- sapply(2:(length(ext.pos) + 1), FUN = function(x){
            return(min(abs(t.ext.density[x] - t.ext.density[x-1]), abs(t.ext.density[x] - t.ext.density[(x+1)])))
        })
        ext <- data.frame(x = ext.pos, y = ext.density, height = ext.height)
        max.ix <- order(ext.density, decreasing = T)
        if(ext.height[max.ix[2]] / ext.height[max.ix[1]] > 0.01){
            cut.df <- ext[c(max.ix[2]:max.ix[1]), ]
            threshold <- x.density$x[cut.df[which.min(cut.df$y), ]$x]
        }else{
            threshold <- NULL
        }
    }else{
        threshold <- NULL
    }

    return(threshold)
}


#
# getBimodalThres <- function(scores){
#     x.density <- density(scores)
#     d.x.density <- diff(x.density$y)
#     d.sign <- (d.x.density > 0) + 0
#
#     ext.pos <- which(d.sign[2:length(d.sign)] - d.sign[1:(length(d.sign)-1)] != 0)
#     if(length(ext.pos) >= 3){
#         ext.density <- x.density$y[ext.pos]
#         t.ext.density <- c(0, ext.density, 0)
#         ext.height <- sapply(2:(length(ext.pos) + 1), FUN = function(x){
#             return(min(abs(t.ext.density[x] - t.ext.density[x-1]), abs(t.ext.density[x] - t.ext.density[(x+1)])))
#         })
#         ext <- data.frame(x = ext.pos, y = ext.density, height = ext.height)
#
#         max.ix <- order(ext.density, decreasing = T)
#         if(ext.height[max.ix[2]] / ext.height[max.ix[1]] > 0.1){
#             cut.df <- ext[c(max.ix[2]:max.ix[1]), ]
#             threshold <- x.density$x[cut.df[which.min(cut.df$y), ]$x]
#         }else{
#             threshold <- NULL
#         }
#     }else{
#         threshold <- NULL
#     }
#     return(threshold)
# }



#' plotMalignancy
#'
#'
#' @return A plot list.
#' @export
#'
plotMalignancy <- function(spot.annotation,
                           coor.names = c("tSNE_1", "tSNE_2"),
                           savePath = NULL){
    ## scatter plot of malignancy
    p.malignType.Point <- pointDRPlot(spot.annotation, value = "Malign.type",
                                      coor.names = coor.names,
                                      colors = c("malignant" = "#f57e87", "nonMalignant" = "#66d5a5"),
                                      legend.position = "right",
                                      legend.title = "Malignancy\n type")

    p.malignScore.Point <- pointDRPlot(spot.annotation, value = "Malign.score",
                                       coor.names = coor.names,
                                       colors = c("white", "#f57e87"),
                                       discrete = F,
                                       limit.quantile = 0.1,
                                       legend.position = "right",
                                       legend.title = "Malignancy\n score")

    p.malignType.bar <- clusterBarPlot(spot.annotation = spot.annotation,
                                       spot.colors = c("malignant" = "#f57e87", "nonMalignant" = "#66d5a5"),
                                       sel.col = "Malign.type",
                                       legend.title = "Malignancy type")

    ## save
    if(!is.null(savePath)){
        ggsave(filename = file.path(savePath, "malignType-point.png"),
               p.malignType.Point, width = 5, height = 3.8, dpi = 500)
        ggsave(filename = file.path(savePath, "malignScore-point.png"),
               p.malignScore.Point, width = 5, height = 3.8, dpi = 500)
        ggsave(filename = file.path(savePath, "malignType-bar.png"),
               p.malignType.bar, width = 6, height = 3, dpi = 500)
    }

    return(list(p.malignType.Point = p.malignType.Point,
                p.malignScore.Point = p.malignScore.Point,
                p.malignType.bar = p.malignType.bar))
}




#' runMalignancy
#'
#' @param expr A Seurat object.
#' @param gene.manifest A data.frame of genes' manifest.
#' @param spot.annotation A data.frame of spots' annotation.
#' @param cutoff The cut-off for min average read counts per gene among
#' reference spots. The default is 0.1.
#' @param minspot An integer number used to filter gene. The default is 3.
#' @param p.value.cutoff The p-value to decide whether the distribution of
#' malignancy score is bimodality.
#' @param ref.data An expression matrix of gene by spot, which is used as the normal reference.
#' The default is NULL, and an immune spots or bone marrow spots expression matrix will be used for human or mouse species, respectively.
#' @param referAdjMat An adjacent matrix for the normal reference data.
#' The larger the value, the closer the spot pair is.
#' The default is NULL, and a SNN matrix of the default ref.data will be used.
#'
#' @return A list of cnvList, reference malignancy score, seurat object,
#' spot.annotatino, bimodal.pvalue, malign.thres, and all generated plots.
#' @export
#'
runMalignancy <- function(SeuratObject,
                          savePath,
                          assay = "Spatial",
                          cutoff = 0.1,
                          minspot = 3,
                          p.value.cutoff = 0.5,
                          coor.names = c("tSNE_1", "tSNE_2"),
                          ref.object = NULL,
                          species = "human",
                          genome = "hg19"){

    gene.manifest <- read.csv(file.path(savePath_data, "geneManifest.txt"),
                              sep = "\t")

    cnvList <- runCNV(SeuratObject,
                      gene.manifest = gene.manifest,
                      assay = assay,
                      cutoff = cutoff, minspot = minspot,
                      ref.object = ref.object,
                      species = species,
                      genome = genome)

    if(is.null(ref.object)){
        if(species == "human"){
            referAdjMat <- readRDS(system.file("rds", "cnvRef_SNN-HM.RDS", package = "stCancer"))
        }else if(species == "mouse"){
            referAdjMat <- readRDS(system.file("rds", "cnvRef_SNN-boneMarrow-MS.RDS", package = "stCancer"))
        }
    }else{
        ref.object <- RunPCA(ref.object, assay = ref.object@active.assay, verbose = FALSE)
        ref.object <- FindVariableFeatures(ref.object, verbose = FALSE)
        ref.object <- FindNeighbors(ref.object, reduction = "pca", dims = 1:30, verbose = FALSE)
        ref.object <- FindClusters(ref.object, resolution = 0.5, verbose = FALSE)
        referAdjMat <- ref.object@graphs$SCT_snn
        colnames(referAdjMat) <- paste0("Normal", 1:ncol(referAdjMat))
        rownames(referAdjMat) <- paste0("Normal", 1:ncol(referAdjMat))
    }
    referScore.smooth <- getMalignScore(cnvList,
                                        "Reference",
                                        method = "smooth",
                                        adjMat = referAdjMat)

    magAdjMat <- SeuratObject@graphs$SCT_snn
    colnames(magAdjMat) <- paste0("Normal", 1:ncol(SeuratObject))
    rownames(magAdjMat) <- paste0("Normal", 1:ncol(SeuratObject))
    obserScore.smooth <- getMalignScore(cnvList,
                                        "Observation",
                                        method = "smooth",
                                        adjMat = magAdjMat)

    up.refer <- quantile(referScore.smooth, 0.995)
    low.refer <- quantile(referScore.smooth, 0.005)
    referScore.smooth <- (referScore.smooth - low.refer) / (up.refer - low.refer)
    obserScore.smooth <- (obserScore.smooth - low.refer) / (up.refer - low.refer)

    all.thres <- getBimodalThres(scores = c(referScore.smooth, obserScore.smooth))
    malign.thres <- getBimodalThres(scores = obserScore.smooth)

    ju.exist.malign <- !is.null(all.thres) | !is.null(malign.thres)

    ## malignancy type
    if(!is.null(all.thres)){
        malign.type <- rep("malignant", length(obserScore.smooth))
        names(malign.type) <- names(obserScore.smooth)
        if(!is.null(malign.thres)){
            malign.type[names(obserScore.smooth)[obserScore.smooth < malign.thres]] <- "nonMalignant"
        }
    }else{
        malign.type <- rep("nonMalignant", length(obserScore.smooth))
        names(malign.type) <- names(obserScore.smooth)
        if(!is.null(malign.thres)){
            malign.type[names(obserScore.smooth)[obserScore.smooth >= malign.thres]] <- "malignant"
        }
    }
    p.malignScore <- malignPlot(obserScore.smooth, referScore.smooth,
                                malign.thres = malign.thres)

    # ## add score and type to spot.annotation
    # spot.annotation$Malign.score <- obserScore.smooth[rownames(spot.annotation)]
    # spot.annotation$Malign.type <- malign.type[rownames(spot.annotation)]
    # expr[["Malign.score"]] <- spot.annotation$Malign.score
    # expr[["Malign.type"]] <- spot.annotation$Malign.type
    #
    # ## plot
    # p.results <- plotMalignancy(spot.annotation = spot.annotation,
    #                             coor.names = coor.names,
    #                             savePath = savePath)
    # p.results[["p.malignScore"]] <- p.malignScore
    # ggsave(filename = file.path(savePath, "malignScore.png"),
    #        p.malignScore, width = 5, height = 4, dpi = 500)

    ## save results
    write.table(cnvList$expr.data[, names(obserScore.smooth)],
                file = file.path(savePath, "inferCNV-observation.txt"),
                quote = F, sep = "\t", row.names = T)
    write.table(cnvList$expr.data[, names(referScore.smooth)],
                file = file.path(savePath, "inferCNV-reference.txt"),
                quote = F, sep = "\t", row.names = T)
    write.table(data.frame(referScore.smooth),
                file = file.path(savePath, "refer-malignScore.txt"),
                quote = F, sep = "\t", row.names = T)
    write.table(data.frame(obserScore.smooth),
                file = file.path(savePath, "obser-malignScore.txt"),
                quote = F, sep = "\t", row.names = T)

    results <- list(
        cnvList = cnvList,
        referScore = referScore.smooth,
        expr = expr,
        ju.exist.malign = ju.exist.malign,
        malign.thres = malign.thres
    )
    return(results)
}


