#' stAnnotation
#'
#' perform statistic
#'
#' @param object A Seurat object
#' @param savePath A path to save the results files(suggest to create a foler named by sample name).
#' @param rm.mito A logical value indicating whether to remove
#' @param crop A logical value indicating whether to crop the image during plotting
#' @param rm.isolated A logical value indicating whether to remove isolated spots
#' @param authorName A character string for authors name and will be shown in the report.
#' @param genReport A logical value indicating whether to generate a .html/.md report (suggest to set TRUE).
#' @param CNV.ref.object A reference seurat object used in inferCNV
#'
#'
#' @return A Seurat object after statistic
#' @export
#'
#' @import Seurat Matrix ggplot2 knitr patchwork
#' @importFrom dplyr "%>%" top_n group_by
#' @importFrom markdown markdownToHTML
#' @importFrom stats median quantile
#' @importFrom utils read.delim read.table write.csv write.table read.csv
stAnnotation <- function(object,
                         savePath,
                         rm.mito = T,
                         rm.ribo = T,
                         # rm.ig = T,
                         assay = NULL,
                         normalization.method = 'SCTransform',
                         topGeneNum = 10,
                         clst.resolution = 0.4,
                         n.markers = 5,
                         rank = 5,
                         species = "human",
                         score.method = 'AddModuleScore',
                         cellTypeScore_fun = "average",
                         cellTypeScore_select = "main",
                         geneSets = NULL,
                         CNV.analysis.func = "copyKAT",
                         CNV.ref.object = NULL,
                         interaction.region.threshold = 20,
                         crop = TRUE,
                         verbose = F,
                         authorName = NULL,
                         genReport = T,
                         bool.NMF = TRUE,
                         bool.CellType = TRUE,
                         bool.CNV = TRUE,
                         bool.interaction = TRUE,
                         bool.tumor.feature = TRUE){

    message("[", Sys.time(), "] START: RUN stAnnotation")

    # collect results and plots
    results.collector <- list()
    plots.collector <- list()

    sampleName <- object@project.name

    savePath <- filePathCheck(savePath)

    savePath <- R.utils::getAbsolutePath(savePath)

    savePath_basic <- file.path(savePath, sampleName)
    savePath_data <- file.path(savePath_basic, "data")
    savePath_cluster <- file.path(savePath_basic, "cluster")
    savePath_pheno <- file.path(savePath_basic, "phenotype")

    if(!dir.exists(savePath_cluster)){
        dir.create(savePath_cluster, recursive = T)
    }

    if(!is.null(assay)){
        assay <- object@active.assay
    }

    # gene.manifest <- read.table(file.path(savePath_data, "geneManifest.txt"),
    #                             header = T, sep = "\t", stringsAsFactors = F)
    # gene.manifest <- subset(gene.manifest, nSpot >= 3)
    # gene.manifest <- subset(gene.manifest, Annotation == "other")

    # remove mito, ribo or ig genes if needed
    if(rm.mito){
        mito.genes <- grep('^MT-', rownames(object), value = TRUE)
        object <- object[!rownames(object) %in% mito.genes, ]
    }
    if(rm.ribo){
        ribo.genes <- grep('^RPL|^RPS|^MRPL|^MRPS', rownames(object), value = TRUE)
        object <- object[!rownames(object) %in% ribo.genes, ]
    }
    # if(rm.ig){
    #   ig.genes <- grep('^IG', rownames(object), value = TRUE)
    #   object <- object[!rownames(object) %in% ig.genes, ]
    # }

    ## ------ Seurat statistic pipeline ------
    if(normalization.method == 'SCTransform'){
        object <- SCTransform(object, assay = object@active.assay, verbose = FALSE)
    }else{
        object <- NormalizeData(object, assay = object@active.assay,
                                normalization.method = normalization.method, verbose = FALSE)
        object <- ScaleData(object, assay = object@active.assay, verbose = FALSE)
        object <- FindVariableFeatures(object, verbose = FALSE)
    }
    object <- RunPCA(object, assay = object@active.assay, verbose = FALSE)
    object <- FindNeighbors(object, reduction = "pca", dims = 1:30, verbose = FALSE)
    object <- FindClusters(object, resolution = clst.resolution, verbose = FALSE)
    suppressWarnings(
        object <- RunUMAP(object, reduction = "pca", dims = 1:30, verbose = F)
    )
    object <- RunTSNE(object, reduction = "pca", dims = 1:30, verbose = F,
                      check_duplicates = FALSE)

    top_gene <- head(VariableFeatures(object), topGeneNum)
    suppressWarnings(
        plot1 <- VariableFeaturePlot(object, cols = c("grey", "#ec7d89"))
    )
    suppressMessages(suppressWarnings(
        hvg <- LabelPoints(plot = plot1, points = top_gene, repel = TRUE) +
            NoLegend()
    ))

    plots.collector[["hvg"]] <- hvg

    ggplot2::ggsave(filename = file.path(savePath_cluster, "hvg.png"),
                    hvg,
                    width = 8,
                    height = 4,
                    dpi = 500)

    clusters <- sort(unique(object@meta.data[['seurat_clusters']]))
    cluster.colors <- getDefaultDimColors(n = length(clusters))
    names(cluster.colors) <- as.character(clusters)

    meta.color <- as.character(object@meta.data[['seurat_clusters']])
    for(i in 1 : length(clusters)){
        meta.color[which(meta.color == clusters[[i]])] <- cluster.colors[[i]]
    }
    object@meta.data[["cluster.color"]] <- meta.color

    p.cls.umap <- Dim_Plot(object,
                           group.by = 'seurat_clusters',
                           reduction = 'umap',
                           legend.name = 'Cluster',
                           legend.color = cluster.colors) +
        ggtitle(NULL)

    p.cls.tsne <- Dim_Plot(object,
                           group.by = 'seurat_clusters',
                           reduction = 'tsne',
                           legend.name = 'Cluster',
                           legend.color = cluster.colors) +
        ggtitle(NULL)

    # suppressWarnings(
    #     p.cls.umap <- pointDRPlot(object,
    #                               feature = 'seurat_clusters',
    #                               reduction.func = "umap",
    #                               colors = cluster.colors,
    #                               base.size = 10,
    #                               legend.title = "Cluster") +
    #         ggplot2::theme(panel.background = element_rect(fill = "white", color = "black"),
    #                        axis.title.x = element_text(size = 20, vjust = -0.2),
    #                        axis.title.y = element_text(size = 20, vjust = 0.2),
    #                        axis.text.x = element_text(size = 16),
    #                        axis.text.y = element_text(size = 16)) +
    #         ggplot2::scale_fill_manual(values = cluster.colors,
    #                                    guide = guide_legend(override.aes = list(size = 8),
    #                                                         keywidth = 0.1,
    #                                                         keyheight = 0.15,
    #                                                         default.unit = "inch"))
    # )

    # suppressWarnings(
    #     p.cls.tsne <- pointDRPlot(object,
    #                               feature = paste0("SCT_snn_res.", clst.resolution),
    #                               reduction.func = "tsne",
    #                               colors = cluster.colors,
    #                               base.size = 10,
    #                               legend.title = "Cluster") +
    #         ggplot2::theme(panel.background = element_rect(fill = "white", color = "black"),
    #                        axis.title.x = element_text(size = 20, vjust = -0.2),
    #                        axis.title.y = element_text(size = 20, vjust = 0.2),
    #                        axis.text.x = element_text(size = 16),
    #                        axis.text.y = element_text(size = 16)) +
    #         ggplot2::scale_fill_manual(values = cluster.colors,
    #                                    guide = guide_legend(override.aes = list(size = 8),
    #                                                         keywidth = 0.1,
    #                                                         keyheight = 0.15,
    #                                                         default.unit = "inch"))
    # )

    p.c <- wrap_plots(p.cls.umap, p.cls.tsne)
    # p.c <- cbind(ggplot2::ggplotGrob(p.cls.umap),
    #              ggplot2::ggplotGrob(p.cls.tsne),
    #              size = "last")
    ggplot2::ggsave(filename = file.path(savePath_cluster, paste0("cluster-dr-", clst.resolution, ".png")),
                    p.c,
                    width = 12,
                    height = 5,
                    dpi = 300)

    plots.collector[["umap"]] <- p.cls.umap
    plots.collector[["tsne"]] <- p.cls.tsne

    p0 <- SpatialFeature_Plot(object,
                              features = 'nCount_Spatial',
                              crop = crop,
                              alpha = c(0,0)) + NoLegend()
    p.cls.sp <- SpatialDim_Plot(object,
                                group.by = 'seurat_clusters',
                                legend.name = 'Cluster',
                                legend.color = cluster.colors,
                                crop = crop)

    # p.cls.sp <- Spatial_Plot(object,
    #                          feature = paste0("SCT_snn_res.", clst.resolution),
    #                          show.tissue = F,
    #                          title = NULL,
    #                          legend.title = "Cluster",
    #                          legend.position = "none",
    #                          colors = cluster.colors,
    #                          discrete = T,
    #                          base.size = 8,
    #                          crop = crop)

    p.sp.c <- wrap_plots(p0, p.cls.sp)
    # p.sp.c <- cbind(ggplot2::ggplotGrob(p0),
    #                 ggplot2::ggplotGrob(p.cls.sp),
    #                 size = "last")
    ggplot2::ggsave(filename = file.path(savePath_cluster, paste0("cluster-sp-", clst.resolution, ".png")),
                    p.sp.c,
                    height = 5,
                    width = 9,
                    dpi = 300)

    plots.collector[["cluster_spatial"]] <- p.cls.sp

    cluster.csv <- data.frame(Barcode = colnames(object),
                              Cluster = object@meta.data[['seurat_clusters']])
    write.csv(cluster.csv,
              file = file.path(savePath_cluster, "cluster.csv"),
              row.names = F,
              quote = F)
    #
    #   saveRDS(object, file.path(savePath_data, "spatial_stat.RDS"))


    ## --------- Differential expression ------------
    diff.expr.genes <- FindAllMarkers(object,
                                      only.pos = TRUE,
                                      min.pct = 0.25,
                                      logfc.threshold = 0.25,
                                      verbose = F)

    diff.expr.genes <- diff.expr.genes[, c("cluster", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
    diff.expr.genes$cluster <- as.numeric(as.character(diff.expr.genes$cluster))
    nCluster <- length(unique(diff.expr.genes$cluster))
    if(!dir.exists(file.path(savePath_cluster, paste0("diff.expr.genes-", clst.resolution)))){
        dir.create(file.path(savePath_cluster, paste0("diff.expr.genes-", clst.resolution)), recursive = T)
    }


    for(ci in 1:nCluster){
        cur.diff.genes <- subset(diff.expr.genes, cluster == unique(diff.expr.genes$cluster)[ci])
        cur.diff.genes <- cur.diff.genes[order(cur.diff.genes$avg_log2FC, decreasing = T), ]
        write.csv(cur.diff.genes,
                  file = file.path(savePath_cluster, paste0("diff.expr.genes-", clst.resolution,
                                                            "/cluster_", unique(diff.expr.genes$cluster)[ci] ,".csv")),
                  quote = F,
                  row.names = F)
    }

    results.collector[["diff.genes"]] <- diff.expr.genes

    top.genes <- diff.expr.genes %>%
        dplyr::group_by(cluster) %>%
        dplyr::top_n(n = n.markers, wt = avg_log2FC)
    top.genes <- top.genes[order(top.genes$cluster, top.genes$avg_log2FC, decreasing = c(F, T)), ]

    de.pre <- preDEheatmap(object = object,
                           col.name = "seurat_clusters",
                           genes = top.genes$gene,
                           slot = "scale.data",
                           min.value = -2.5,
                           max.value = 2.5)

    # ann_colors <- unique(object@meta.data[["cluster.color"]])
    # ann_colors <- ann_colors[order(as.numeric(unique(object@meta.data[["seurat_clusters"]])))]
    # names(ann_colors) <- sort(unique(object@meta.data[["seurat_clusters"]]))
    ann_colors <- list(cluster.colors)
    names(ann_colors) <- "clusters"
    names(de.pre$spot.cluster) <- "clusters"

    p.de.heatmap <- pheatmap::pheatmap(de.pre$expr.data,
                                       color = colorRampPalette(c("#4393C3", "white", "#D6604D"))(100),
                                       annotation_col = de.pre$spot.cluster,
                                       annotation_colors = ann_colors,
                                       fontsize = 7,
                                       gaps_col = de.pre$gaps_col,
                                       cluster_rows = F,
                                       cluster_cols = F,
                                       show_colnames = F,
                                       use_raster = F)

    DEplot.height <- 0.5 + 0.1 * n.markers * length(unique(object@meta.data[["seurat_clusters"]]))
    ggplot2::ggsave(filename = file.path(savePath_cluster, paste0("DE-heatmap-", clst.resolution, ".png")),
                    p.de.heatmap,
                    width = 8,
                    height = DEplot.height,
                    dpi = 500)

    plots.collector[["diff.genes"]] <- p.de.heatmap

    ## --------- NMF analysis ------------
    if(bool.NMF){
        savePath_NMF <- file.path(savePath_basic, "NMF")

        if(!dir.exists(savePath_NMF)){
            dir.create(savePath_NMF, recursive = T)
        }

        results.NMF <- runNMF(object = object,
                              savePath = savePath_NMF,
                              rank = rank)
        object <- results.NMF$object
        results.collector[["NMF"]] <- results.NMF$results
        plots.collector[["NMF"]] <- results.NMF$plots
        bool.NMF <- results.NMF$bool.completed
    }

    ## --------- Cell type score ------------
    if(bool.CellType){
        savePath_Score <- file.path(savePath_basic, "cell_type")

        if(!dir.exists(savePath_Score)){
            dir.create(savePath_Score, recursive = T)
        }

        results.CellType <- cellTypeScore(object,
                                          savePath = savePath_Score,
                                          geneSets = geneSets,
                                          method = cellTypeScore_fun,
                                          select = cellTypeScore_select,
                                          species = species,
                                          crop = crop)
        object <- results.CellType$object
        results.collector[["cell_type"]] <- results.CellType$results
        plots.collector[["cell_type"]] <- results.CellType$plots
        bool.CellType <- results.CellType$bool.completed
    }


    ## --------- CNV analysis ------------
    if(bool.CNV){
        savePath_malig <- file.path(savePath_basic, "malignancy")

        if(!dir.exists(savePath_malig)){
            dir.create(savePath_malig, recursive = T)
        }

        results.CNV <- runCNVAnalysis(object,
                                      savePath_malig,
                                      species = species,
                                      analysis.func = CNV.analysis.func,
                                      ref.object = CNV.ref.object)
        results.collector[["CNV"]] <- results.CNV$results
        plots.collector[["CNV"]] <- results.CNV$plots
        bool.CNV <- results.CNV$bool.completed
    }


    ## --------- Interaction ------------
    if(bool.interaction){
        savePath_interact <- file.path(savePath_basic, "interact")

        if(!dir.exists(savePath_interact)){
            dir.create(savePath_interact, recursive = T)
        }

        results.interaction <- SpatialInteraction(object =  object,
                                                  savePath = savePath_interact,
                                                  # genePath = file.path(savePath_data, "filtered_feature_bc_matrix/features.tsv.gz"),
                                                  species = species,
                                                  region.threshold = interaction.region.threshold)
        results.collector[["interaction"]] <- results.interaction$results
        plots.collector[["interaction"]] <- results.interaction$plots
        bool.interaction <- results.interaction$bool.completed
    }

    # ## --------- EMT ---------
    # if(bool.EMT){
    #     savePath_pheno <- file.path(savePath_basic, "phenotype")
    #
    #     if(!dir.exists(savePath_pheno)){
    #         dir.create(savePath_pheno, recursive = T)
    #     }
    #
    #     results.EMT <- EMTAnalysis(object,
    #                                savePath = savePath_pheno,
    #                                species = species)
    #     object <- results.EMT$object
    #     results.collector[["EMT"]] <- results.EMT$results
    #     plots.collector[["EMT"]] <- results.EMT$plots
    #     bool.EMT <- results.EMT$bool.completed
    # }
    #
    # ## --------- Cell cycle ---------
    # if(bool.CellCycle){
    #     savePath_pheno <- file.path(savePath_basic, "phenotype")
    #
    #     if(!dir.exists(savePath_pheno)){
    #         dir.create(savePath_pheno, recursive = T)
    #     }
    #     results.CellCycle <- CellCycleAnalysis(object,
    #                                            savePath = savePath_pheno,
    #                                            species = species)
    #     object <- results.CellCycle$object
    #     results.collector[["CellCycle"]] <- results.CellCycle$results
    #     plots.collector[["CellCycle"]] <- results.CellCycle$plots
    #     bool.CellCycle <- results.CellCycle$bool.completed
    # }
    #
    # ## --------- Stem ---------
    # if(bool.stem){
    #     savePath_pheno <- file.path(savePath_basic, "phenotype")
    #
    #     if(!dir.exists(savePath_pheno)){
    #         dir.create(savePath_pheno, recursive = T)
    #     }
    #     results.stem <- StemnessAnalysis(object,
    #                                      savePath = savePath_pheno,
    #                                      species = species)
    #     object <- results.stem$object
    #     results.collector[["stem"]] <- results.stem$results
    #     plots.collector[["stem"]] <- results.stem$plots
    #     bool.stem <- results.stem$bool.completed
    # }

    if(bool.tumor.feature){
        ## --------- CancerSEA characters --------
        savePath_cancerSEA <- file.path(savePath_basic, 'cancer_state')

        if(!dir.exists(savePath_cancerSEA)){
            dir.create(savePath_cancerSEA, recursive = T)
        }
        results.cancer <- TumorCharacters(object,
                                          savePath = savePath_cancerSEA,
                                          species = species,
                                          method = score.method,
                                          crop = crop)
        object <- results.cancer$object
        results.collector[["Tumor"]] <- results.cancer$results
        plots.collector[["Tumor"]] <- results.cancer$plots

        ## --------- TLS analysis --------
        savePath_cancerSEA <- file.path(savePath_basic, 'cancer_state')

        if(!dir.exists(savePath_cancerSEA)){
            dir.create(savePath_cancerSEA, recursive = T)
        }
        results.cancer <- TLSAnalysis(object,
                                      savePath = savePath_cancerSEA,
                                      species = species,
                                      method = score.method,
                                      crop = crop)
        object <- results.cancer$object
        results.collector[["TLS"]] <- results.cancer$results
        plots.collector[["TLS"]] <- results.cancer$plots
    }

    ## --------- Save data ---------
    saveRDS(object, file.path(savePath_data, "spatial_anno.RDS"))

    ## --------- HTML report ------------
    if(genReport){
        knitr::knit(file.path(system.file(package = "stCancer"), "rmd/stAnno.Rmd"),
                    file.path(savePath_basic, "report-stAnno.md"))
        markdown::markdownToHTML(file.path(savePath_basic, 'report-stAnno.md'),
                                 file.path(savePath_basic, 'report-stAnno.html'))
    }

    return(list(object = object,
                results = results.collector,
                plots = plots.collector))
}
