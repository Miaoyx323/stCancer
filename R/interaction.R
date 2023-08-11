#' regionprops
#'
#' find all the connected domains of `cluster`
#'
#' @param object A Seurat object
#' @param cluster The cluster to analyze
#' @param cluster.info The label to store cluster information in `object@meta.data`
#' @param region.threshold Only domains whose area are larger than `region.threshold` are remained.
#' @param slice
#'
#'
#' @return A list of connected domains with information of barcodes and clusters
#' @export
#'
#' @importFrom stats na.omit
regionprops <- function(object,
                        cluster,
                        region.select.method = "threshold",
                        region.threshold = 20,
                        cluster.info = "seurat_clusters",
                        slice = "slice1"){
    # get the data.frame that annotates the cluster of each spot
    coordinates <- object@images[[slice]]@coordinates
    coordinates <- cbind(coordinates, as.data.frame(object[[cluster.info]]))

    basic_coordinates <- coordinates[colnames(object[, object[[cluster.info]] == cluster]), ]

    clusters <- unique(object@meta.data[[cluster.info]])
    clusters <- as.numeric(as.character(clusters))

    # record the results
    results.data.frame <- list()
    i.results.data.frame <- 1

    for(i in 1 : length(clusters)){
        i_cluster <- clusters[[i]]
        if(i_cluster == cluster){
            next
        }
        # find all the area of cluster and i_cluster
        i_coordinates <- coordinates[colnames(object[, object[[cluster.info]] == i_cluster]), ]
        i_coordinates <- rbind(i_coordinates, basic_coordinates)
        # find all the boundary spots
        i_coordinates$label <- -1
        for(j in 1 : nrow(i_coordinates)){
            j_row <- i_coordinates[j, ]$row
            j_col <- i_coordinates[j, ]$col
            j_barcode <- rownames(i_coordinates[j, ])
            j_cluster <- as.numeric(as.character(i_coordinates[j, cluster.info]))
            nei <- get_neighbors_indices(j_row, j_col)
            for(k in 1 : length(nei)){
                k_coordinates <- i_coordinates[i_coordinates$row == nei[[k]][1] &
                                                   i_coordinates$col == nei[[k]][2], ]
                if(nrow(na.omit(k_coordinates)) == 0){
                    next
                }
                # if the cluster of its neighbors is different, it is on the boundary
                if(as.numeric(as.character(k_coordinates[1, cluster.info])) != j_cluster){
                    i_coordinates[j_barcode, "label"] <- 1
                }
            }
        }
        # get all the boundary
        i_boundary <- i_coordinates %>% subset(label == 1)
        # no connection
        if(nrow(na.omit(i_boundary)) == 0){
            next
        }
        # Connected domain labeling
        i_boundary$label <- -1
        label_id <- 1
        for(j in 1 : nrow(i_boundary)){
            if(i_boundary[j, "label"] != -1){
                next
            }

            i_boundary[j, "label"] <- label_id

            j_row <- i_boundary[j, ]$row
            j_col <- i_boundary[j, ]$col
            j_barcode <- rownames(i_boundary[j, ])
            nei <- get_neighbors_indices(j_row, j_col)
            # area list
            spot_list <- list()
            spot_list[[1]] <- c(j_row, j_col)
            while(length(spot_list) != 0){
                ii_spot <- spot_list[[1]]
                spot_list[[1]] <- NULL
                nei <- get_neighbors_indices(ii_spot[1], ii_spot[2])
                for(k in 1 : length(nei)){
                    k_boundary <- i_boundary[i_boundary$row == nei[[k]][1] &
                                                 i_boundary$col == nei[[k]][2], ]
                    if(nrow(na.omit(k_boundary)) == 0){
                        next
                    }
                    if(k_boundary$label == -1){
                        spot_list[[length(spot_list) + 1]] <- c(k_boundary$row, k_boundary$col)
                        i_boundary[i_boundary$row == nei[[k]][1] &
                                       i_boundary$col == nei[[k]][2], "label"] <- label_id
                    }
                }
            }
            label_id <- label_id + 1
        }
        # judge the size of area
        for(j in 1 : label_id){
            j_boundary <- i_boundary %>% subset(label == j)
            j_boundary <- na.omit(j_boundary)
            if(nrow(j_boundary) >= region.threshold & nrow(j_boundary) != 0){
                # add to results.data.frame
                results.data.frame[[i.results.data.frame]] <- data.frame(barcode = rownames(j_boundary),
                                                                         group = j_boundary[[cluster.info]])
                i.results.data.frame <- i.results.data.frame + 1
            }
        }
    }

    return(results.data.frame)

}






#' SpatialInteraction
#'
#' interaction
#'
#' @param object A Seurat object
#' @param savePath A path to store results
#' @param genePath The path of features.tsv.gz
#' @param analysis.type "biggest" means the function will analyze the interactions of the biggest cluster.
#' "single" means the function will analyze the interaction of cluster `cluster`.
#' "all" means the function will analyze all the clusters
#' @param cluster Only when `analysis.type` is "single", this parameter means the cluster to analyze
#'
#'
#' @import org.Hs.eg.db org.Mm.eg.db
#' @import clusterProfiler
#' @return
#' @export
#'
SpatialInteraction <- function(object,
                               savePath,
                               # genePath, # path of features.tsv.gz
                               analysis.type = "biggest", # biggest, single, all
                               cluster = NULL,
                               area.data.frame = NULL, # area barcodes and their labels
                               species = "human",
                               label_key = 'seurat_clusters',
                               ...){

    # collect results and plots
    results.collector <- list()
    plots.collector <- list()

    sampleName <- object@project.name

    if(species == 'human'){
        OrgDb.use <- org.Hs.eg.db
    }else if(species == 'mouse'){
        OrgDb.use <- org.Mm.eg.db
    }else{
        stop(paste0("The '", species, "' is not supported in the interaction part."))
        return(list(results = results.collector,
                    plots = plots.collector,
                    bool.completed = FALSE))
    }

    genes <- rownames(object)
    genes <- bitr(genes, fromType = 'SYMBOL',
                  toType = 'ENSEMBL',
                  OrgDb = OrgDb.use)

    # genes <- read.table(gzfile(genePath, "rt"))
    # colnames(genes) <- c("A", "B", "C", "D")

    if(analysis.type == "biggest"){
        clusters.info <- object@meta.data[[label_key]]
        clusters <- unique(clusters.info)
        clusters <- as.character(clusters)
        # clusters <- as.numeric(as.character(clusters))
        n_each_cluster <- matrix(nrow = 1,
                                 ncol = length(clusters))
        for(i in 1 : length(clusters)){
            n_each_cluster[i] <- sum(clusters.info == clusters[i])
        }
        cluster <- clusters[which(n_each_cluster == max(n_each_cluster))]
        prepareCellPhoneDB(object,
                           cluster,
                           savePath,
                           genes,
                           ...)
    }else if(analysis.type == "all"){
        # get all the clusters
        clusters <- unique(object@meta.data[[label_key]])
        clusters <- as.numeric(as.character(clusters))
        for(i in 1 : length(clusters)){
            cluster <- clusters[i]
            prepareCellPhoneDB(object,
                               cluster,
                               savePath,
                               genes,
                               ...)
        }
    }else if(analysis.type == "single"){
        if(is.null(cluster)){
            clusters.info <- object@meta.data[[label_key]]
            clusters <- unique(clusters.info)
            clusters <- as.numeric(as.character(clusters))
            n_each_cluster <- matrix(nrow = 1,
                                     ncol = length(clusters))
            for(i in 1 : length(clusters)){
                n_each_cluster[i] <- sum(clusters.info == clusters[i])
            }
            cluster <- clusters[which(n_each_cluster == max(n_each_cluster))]
            prepareCellPhoneDB(object,
                               cluster,
                               savePath,
                               genes,
                               ...)
        }else{
            prepareCellPhoneDB(object,
                               cluster,
                               savePath,
                               genes,
                               ...)
        }
    }else{

    }

    results.collector <- callCellPhoneDB(object,
                                         savePath)

    labels <- list.files(savePath)
    files <- file.path(savePath, labels)

    for(i in 1 : length(files)){
        outPath <- file.path(files[[i]], "out")
        if(dir.exists(outPath)){
            p <- InteractionDotPlot(outPath)
            p0 <- png::readPNG(source = file.path(files[[i]], 'area.png'))
            plots.collector[[paste0("area ", i)]] <- list(area = p0,
                                                          interaction = p)

            knitr::knit(file.path(system.file(package = "stCancer"), "rmd/interaction.Rmd"),
                        file.path(files[[i]], "interaction.md"))
            markdown::markdownToHTML(file.path(files[[i]], 'interaction.md'),
                                     file.path(files[[i]], 'interaction.html'))
        }
    }

    return(list(results = results.collector,
                plots = plots.collector,
                bool.completed = TRUE))
}




InteractionDotPlot <- function(outPath,
                               n.top = NULL,
                               only.inter = F){

    signif_level <- c(`***` = 0.001, `**` = 0.01, `*` = 0.05, ns = 1)

    # read the output files of CellPhoneDB
    mean.df <- read.table(file.path(outPath, "means.txt"),
                          header = T,
                          sep = "\t",
                          check.names = F)
    rownames(mean.df) <- mean.df$interacting_pair

    pvalue.df <- read.table(file.path(outPath, "pvalues.txt"),
                            header = T,
                            sep = "\t",
                            check.names = F)
    rownames(pvalue.df) <- pvalue.df$interacting_pair

    if(only.inter){
        mean.df <- mean.df[, c(1:11, 13, 14)]
        pvalue.df <- pvalue.df[, c(1:11, 13, 14)]
    }

    if(is.null(n.top)){
        # select the pairs whose interaction is significant, p-value < 0.05
        pvalue.df$Exist_signif <- rowSums(pvalue.df[, c(12:ncol(pvalue.df))] < 0.05)
        pvalue.df <- pvalue.df[pvalue.df$Exist_signif > 0, ]
        pvalue.df <- pvalue.df[, 1:(ncol(pvalue.df) - 1)]
    }else{
        row.mins <- vector()
        for(i in 1 : nrow(pvalue.df)){
            row.mins[i] <- min(pvalue.df[i, c(12:ncol(pvalue.df))])
        }
        idx <- order(row.mins)[1:n.top]
        pvalue.df <- pvalue.df[idx, ]
    }


    mean.df <- mean.df[rownames(pvalue.df), ]
    if(!all.equal(rownames(mean.df), rownames(pvalue.df))){
        stop("Not equal!")
    }

    me.df <- mean.df[, 12:ncol(mean.df)]
    colnames(me.df) <- paste0(colnames(me.df), "_Mean")
    pv.df <- pvalue.df[, 12:ncol(pvalue.df)]
    colnames(pv.df) <- paste0(colnames(pv.df), "_pvalue")
    col.order <- order(c(1:ncol(me.df)*2-1, 1:ncol(me.df)*2))
    info.df <- cbind(me.df, pv.df)
    info.df <- info.df[, col.order]
    info.df <- cbind(mean.df[, 1:11], info.df)

    mean.melt <- melt(mean.df[, c(2, 12:ncol(mean.df))])
    pvalue.melt <- melt(pvalue.df[, c(2, 12:ncol(mean.df))])
    check.equal <- all.equal(paste(mean.melt$interacting_pair, mean.melt$variable),
                             paste(pvalue.melt$interacting_pair, pvalue.melt$variable))

    if(check.equal){
        p.df <- cbind(mean.melt, pvalue.melt$value)
        colnames(p.df) <- c("interacting_pair", "cluster_pair", "Mean", "p.value")
        p.df$logMean <- log2(p.df$Mean + 1)
        p.df$log.pvalue <- -log10(p.df$p.value)
        p.df$signif_level <- sapply(1:nrow(p.df), function(i){
            names(signif_level)[min(which(p.df$p.value[i] <= signif_level))]
        })
        p.df$signif_level <- factor(p.df$signif_level, levels = c("ns", "*", "**", "***"))

        p <- ggplot(p.df, aes(x = interacting_pair,
                              y = cluster_pair,
                              color = logMean,
                              size = signif_level)) +
            geom_point() +
            labs(x = NULL, y = NULL, color = "Mean", size = "Significance") +
            scale_color_gradientn(colours = rev(brewer.pal(11, "Spectral"))) +
            scale_size_manual(values = c("ns" = 1, "*" = 3, "**" = 4, "***" = 5)) +
            theme_bw() +
            theme(axis.text.y = element_text(angle = 45, vjust = 0.9, hjust = 0.9,
                                             color = 'black'),
                  axis.text.x = element_text(angle = -90, color = 'black'))

        ggsave(file.path(outPath, "dot.png"), limitsize = FALSE,
               p, height = 5, width = 0.5+nrow(pvalue.df)/6.0, dpi = 300)
    }
    return(p)
}

prepareCellPhoneDB <- function(object,
                               cluster,
                               savePath,
                               genes,
                               label_key = 'seurat_clusters',
                               ...){

    results <- regionprops(object,
                           cluster = cluster,
                           cluster.info = label_key,
                           ...)

    if(is.null(results) | length(results) == 0){
        return(0)
    }

    for(i in 1 : length(results)){
        savePath_cellphonedb <- file.path(savePath, paste0("cluster_", cluster, "_region_", i))
        if(!dir.exists(savePath_cellphonedb)){
            dir.create(savePath_cellphonedb, recursive = T)
        }

        barcodes <- results[[i]]
        colnames(barcodes) <- c("barcode", "group")
        write.table(barcodes,
                    file.path(savePath_cellphonedb,
                              "meta.txt"),
                    row.names = F,
                    quote = F,
                    sep = "\t")

        expr.data <- GetAssayData(object, slot = "counts")[, barcodes$barcode]
        expr.data <- apply(expr.data, 2, function(x) (x/sum(x))*10000)

        # gene <- as.data.frame(rownames(expr.data))
        # colnames(gene) <- "Gene"

        genes <- genes[genes$SYMBOL %in% rownames(expr.data), ]
        genes <- genes[!duplicated(genes$ENSEMBL), ]
        expr.data <- expr.data[genes$SYMBOL, ]
        rownames(expr.data) <- genes$ENSEMBL
        expr.data <- as.data.frame(expr.data)

        # com.genes <- intersect(gene$Gene, genes$B)
        # tmpp <- genes %>% subset(B %in% com.genes)
        # tmppp <- tmpp[!duplicated(tmpp$B), ]
        #
        # expr.data <- cbind.data.frame(gene, expr.data)
        # expr.data <- expr.data %>% subset(Gene %in% com.genes)
        # expr.data$Gene <- tmppp$A
        #
        # rownames(expr.data) <- expr.data$Gene
        # expr.data <- expr.data %>% subset(select = -c(Gene))

        write.table(expr.data,
                    file.path(savePath_cellphonedb,
                              "counts.txt"),
                    row.names = T,
                    quote = F,
                    sep = "\t")

        # draw the interaction area
        object$barcode <- colnames(object)
        subobject <- subset(object, barcode %in% barcodes$barcode)

        # subobject <- object[, barcode$barcode]
        ggplot2::ggsave(filename = file.path(savePath_cellphonedb, "area.png"),
                        SpatialDim_Plot(subobject,
                                        label_key,
                                        crop = F,
                                        legend.color =
                                            getDefaultDimColors(n = 2)) +
                            theme(legend.position = 'top'),
                        # Spatial_Plot(subobject,
                        #              "seurat_clusters",
                        #              show.tissue = T,
                        #              crop = T,
                        #              colors = "cluster.color",
                        #              discrete = T,
                        #              pt.size = 1.6,
                        #              base.size = 8),
                        width = 6,
                        height = 5,
                        dpi = 300)
    }
}




callCellPhoneDB <- function(object,
                            savePath,
                            pvalues.threshold = 0.01){

    # collect results
    results.collector <- list()

    sampleName <- object@project.name

    labels <- list.files(savePath)
    files <- file.path(savePath, labels)

    cellphone.db <- readRDS(file.path(system.file(package = "stCancer"), "rds/cellphonedb.RDS"))

    for(i in 1 : length(files)){
        meta <- read.table(file.path(files[[i]], "meta.txt"),
                           sep = "\t",
                           header = T)
        counts <- read.table(file.path(files[[i]], "counts.txt"),
                             sep = "\t",
                             header = T,
                             row.names = 1)
        results <- myCellPhoneDB(meta,
                                 counts,
                                 files[[i]],
                                 cellphone.db,
                                 pvalues.threshold)

        results.collector[[paste0("area ", i)]] <- list(meta = meta,
                                                        counts = counts,
                                                        means = results$means,
                                                        pvalues = results$pvalues)

    }

    return(results.collector)
}



myCellPhoneDB <- function(meta,
                          counts,
                          savePath,
                          cellphone.db,
                          pvalues.threshold = 0.01){

    ### transform mRNA counts to protein counts
    # read all tables
    complex_table <- cellphone.db[["complex_table"]]
    complex_composition_table <- cellphone.db[["complex_composition_table"]]
    gene_table <- cellphone.db[["gene_table"]]
    interaction_table <- cellphone.db[["interaction_table"]]
    multidata_table <- cellphone.db[["multidata_table"]]
    protein_table <- cellphone.db[["protein_table"]]

    # find common genes
    common_gene <- intersect(rownames(counts), gene_table$ensembl)

    # filter data and tables
    counts <- counts[common_gene, ]
    gene_table <- gene_table %>%
        subset(ensembl %in% common_gene)

    ## get the gene and protein table, and find out all the possible protein
    gene_protein_table <- merge(gene_table,
                                protein_table,
                                by.x = "protein_id",
                                by.y = "id_protein")
    # find single protein
    protein_multidata_id <- gene_protein_table$protein_multidata_id

    # find complex
    filtered_complex_composition_table <- merge(complex_composition_table,
                                                gene_protein_table,
                                                by.x = "protein_multidata_id",
                                                by.y = "protein_multidata_id")
    complex_multidata_id <- unique(filtered_complex_composition_table$complex_multidata_id)
    for(i in 1 : length(complex_multidata_id)){
        i_id <- complex_multidata_id[i]
        i_table <- filtered_complex_composition_table %>%
            subset(complex_multidata_id == i_id)
        if(nrow(i_table) != unique(i_table$total_protein)){
            complex_multidata_id[i] <- NA
        }
    }
    complex_multidata_id <- na.omit(complex_multidata_id)

    # combine single protein and complex
    multidata_id <- c(protein_multidata_id, complex_multidata_id)

    ## build protein counts
    protein_counts <- matrix(nrow = length(multidata_id),
                             ncol = ncol(counts))
    rownames(protein_counts) <- as.character(multidata_id)
    colnames(protein_counts) <- colnames(counts)

    gene_protein_table <- gene_protein_table[!duplicated(gene_protein_table$protein_multidata_id), ]

    rownames(gene_protein_table) <- gene_protein_table$protein_multidata_id
    gene_protein_table <- gene_protein_table[as.character(protein_multidata_id), ]

    for(i in 1 : length(protein_multidata_id)){
        i_ensembl <- gene_protein_table$ensembl[i]
        protein_counts[i, ] <- as.matrix(counts[i_ensembl, ])
    }

    ## build complex counts
    for(i in 1 : length(complex_multidata_id)){
        i_table <- complex_composition_table[
            complex_composition_table$complex_multidata_id %in%
                complex_multidata_id[i], ]
        i_protein_ids <- i_table$protein_multidata_id
        i_gene_protein_table <- gene_protein_table[as.character(i_protein_ids), ]
        i_ensembl <- i_gene_protein_table$ensembl
        protein_counts[as.character(complex_multidata_id[i]), ] <-
            cal_complex(counts[i_ensembl, ])
        protein_counts[as.character(complex_multidata_id[i]), ] <-
            cal_complex(counts[i_ensembl, ])
    }

    # find possible interaction
    filtered_interaction_table <- interaction_table %>%
        subset(multidata_1_id %in% multidata_id) %>%
        subset(multidata_2_id %in% multidata_id)

    ## statistic method
    statistic_results <- statistic_method(meta,
                                          protein_counts,
                                          filtered_interaction_table)

    ## get pvalue
    clusters <- unique(meta[, 2])
    cluster_1 <- clusters[1]
    cluster_2 <- clusters[2]
    n_interaction <- nrow(filtered_interaction_table)

    protein_counts_1 <- protein_counts[, make.names(meta[meta[, 2] == cluster_1, ][, 1])]
    protein_counts_2 <- protein_counts[, make.names(meta[meta[, 2] == cluster_2, ][, 1])]
    mean_protein_counts_1 <- rowMeans(protein_counts_1)
    mean_protein_counts_2 <- rowMeans(protein_counts_2)

    # build real matrix
    real_mean <- matrix(nrow = n_interaction,
                        ncol = 4)
    colnames(real_mean) <- c(paste0(cluster_1, "|", cluster_1),
                             paste0(cluster_1, "|", cluster_2),
                             paste0(cluster_2, "|", cluster_1),
                             paste0(cluster_2, "|", cluster_2))
    rownames(real_mean) <- filtered_interaction_table$id_cp_interaction


    for(i in 1 : n_interaction){
        real_mean[i, paste0(cluster_1, "|", cluster_1)] <-
            cal_mean(mean_protein_counts_1[as.character(filtered_interaction_table[i, ]$multidata_1_id)],
                     mean_protein_counts_1[as.character(filtered_interaction_table[i, ]$multidata_2_id)])
        real_mean[i, paste0(cluster_1, "|", cluster_2)] <-
            cal_mean(mean_protein_counts_1[as.character(filtered_interaction_table[i, ]$multidata_1_id)],
                     mean_protein_counts_2[as.character(filtered_interaction_table[i, ]$multidata_2_id)])
        real_mean[i, paste0(cluster_2, "|", cluster_1)] <-
            cal_mean(mean_protein_counts_2[as.character(filtered_interaction_table[i, ]$multidata_1_id)],
                     mean_protein_counts_1[as.character(filtered_interaction_table[i, ]$multidata_2_id)])
        real_mean[i, paste0(cluster_2, "|", cluster_2)] <-
            cal_mean(mean_protein_counts_2[as.character(filtered_interaction_table[i, ]$multidata_1_id)],
                     mean_protein_counts_2[as.character(filtered_interaction_table[i, ]$multidata_2_id)])
    }

    # pvalue matrix
    pvalue_matrix <- matrix(nrow = n_interaction,
                            ncol = 4)
    colnames(pvalue_matrix) <- c(paste0(cluster_1, "|", cluster_1),
                                 paste0(cluster_1, "|", cluster_2),
                                 paste0(cluster_2, "|", cluster_1),
                                 paste0(cluster_2, "|", cluster_2))
    rownames(pvalue_matrix) <- filtered_interaction_table$id_cp_interaction

    for(i in 1 : 4){
        for(j in 1 : n_interaction){
            n_bigger <- 0
            for(iter in 1 : 1000){
                if(statistic_results[[i]][j, iter] > real_mean[j, i]){
                    n_bigger <- n_bigger + 1
                }
            }
            pvalue_matrix[j, i] <- n_bigger / 1000
        }
    }

    real_mean <- real_mean[which(rowSums(real_mean) > 0), ]
    pvalue_matrix <- pvalue_matrix[rownames(real_mean), ]

    # select interactions whose pvalue is less than pvalues.threshold
    pvalue_matrix <- pvalue_matrix[rowSums(pvalue_matrix <= pvalues.threshold) != 0, ]
    valid_interaction <- rownames(pvalue_matrix)

    if(is.null(valid_interaction)){
        return(1)
    }

    mean_matrix <- real_mean[valid_interaction, ]
    pvalue_matrix <- pvalue_matrix[valid_interaction, ]

    rownames(filtered_interaction_table) <- filtered_interaction_table$id_cp_interaction
    filtered_interaction_table <- filtered_interaction_table[valid_interaction, ]

    ## get table base
    table_base <- build_table_base(filtered_interaction_table = filtered_interaction_table,
                                   multidata_table = multidata_table,
                                   gene_protein_table = gene_protein_table)

    ## means.txt
    means.txt <- cbind(table_base,
                       as.data.frame(mean_matrix))
    pvalues.txt <- cbind(table_base,
                         as.data.frame(pvalue_matrix))

    savePath <- file.path(savePath, "out")
    dir.create(savePath, recursive = T)

    write.table(means.txt,
                file.path(savePath,
                          "means.txt"),
                row.names = F,
                quote = F,
                sep = "\t")

    write.table(pvalues.txt,
                file.path(savePath,
                          "pvalues.txt"),
                row.names = F,
                quote = F,
                sep = "\t")

    return(list(means = means.txt,
                pvalues = pvalues.txt))
}


statistic_method <- function(meta,
                             counts,
                             interaction_table,
                             iteration = 1000){

    clusters <- unique(meta[, 2])

    if(length(clusters) != 2){
        stop("The number of types of meta is not 2!")
    }

    cluster_1 <- clusters[1]
    cluster_2 <- clusters[2]

    n_cluster <- nrow(meta)
    n_cluster_1 <- sum(meta[, 2] == cluster_1)
    n_cluster_2 <- sum(meta[, 2] == cluster_2)
    n_interaction <- nrow(interaction_table)

    # "results" includes 4 matrices, and the name of matrix stands for "ligand_receptor"
    results <- list()
    results[[paste0(cluster_1, "|", cluster_1)]] <- matrix(nrow = n_interaction,
                                                           ncol = iteration)
    results[[paste0(cluster_1, "|", cluster_2)]] <- matrix(nrow = n_interaction,
                                                           ncol = iteration)
    results[[paste0(cluster_2, "|", cluster_1)]] <- matrix(nrow = n_interaction,
                                                           ncol = iteration)
    results[[paste0(cluster_2, "|", cluster_2)]] <- matrix(nrow = n_interaction,
                                                           ncol = iteration)

    for(i in 1 : iteration){
        rank <- sample(c(1:n_cluster), n_cluster)
        i_counts_cluster_1 <- counts[, rank[1:n_cluster_1]]
        i_counts_cluster_2 <- counts[, rank[(n_cluster_1+1):n_cluster]]

        for(j in 1 : nrow(interaction_table)){
            j_table <- interaction_table[j, ]
            interaction_1_l <- mean(i_counts_cluster_1[as.character(j_table$multidata_1_id), ])
            interaction_1_r <- mean(i_counts_cluster_1[as.character(j_table$multidata_2_id), ])
            interaction_2_l <- mean(i_counts_cluster_2[as.character(j_table$multidata_1_id), ])
            interaction_2_r <- mean(i_counts_cluster_2[as.character(j_table$multidata_2_id), ])

            results[[paste0(cluster_1, "|", cluster_1)]][j, i] <- cal_mean(interaction_1_l, interaction_1_r)
            results[[paste0(cluster_1, "|", cluster_2)]][j, i] <- cal_mean(interaction_1_l, interaction_2_r)
            results[[paste0(cluster_2, "|", cluster_1)]][j, i] <- cal_mean(interaction_2_l, interaction_1_r)
            results[[paste0(cluster_2, "|", cluster_2)]][j, i] <- cal_mean(interaction_2_l, interaction_2_r)
        }
    }

    return(results)
}

cal_mean <- function(mean_1,
                     mean_2){
    if(mean_1 == 0 | mean_2 == 0){
        return(0)
    }else{
        return((mean_1 + mean_2) / 2)
    }
}

cal_complex <- function(counts){
    if(sum(counts == 0) != 0){
        return(0)
    }else{
        # return(colMeans(counts))
        return(apply(counts, 2, min))
    }
}


build_table_base <- function(filtered_interaction_table,
                             multidata_table,
                             gene_protein_table){
    ## id_cp_interaction
    id_cp_interaction <- filtered_interaction_table$id_cp_interaction
    ## interaction_pair, ensembl, partner, secreted and is_integrin
    ligand_id <- filtered_interaction_table$multidata_1_id
    receptor_id <- filtered_interaction_table$multidata_2_id
    interaction_pair <- list()
    partner_a <- list()
    partner_b <- list()
    ensembl_a <- list()
    ensembl_b <- list()
    secreted <- list()
    receptor_a <- list()
    receptor_b <- list()
    is_integrin <- list()

    for(i in 1 : length(ligand_id)){
        # secreted flag
        flag_secreted <- FALSE
        # is_integrin flag
        flag_integrin <- FALSE
        # ligand
        i_table <- subset(gene_protein_table, protein_multidata_id == as.character(ligand_id[i]))
        # complex
        if(nrow(i_table) != 1){
            i_table <- subset(multidata_table, id_multidata == as.character(ligand_id[i]))
            i_ligand <- i_table$name
            partner_a[[i]] <- paste0("complex:", i_table$name)
            ensembl_a[[i]] <- ""
            receptor_a[[i]] <- i_table$receptor == 1

            if(i_table$secreted){
                flag_secreted <- TRUE
            }
            flag_integrin <- TRUE
        }else{ # single protein
            i_ligand <- i_table$gene_name
            ensembl_a[[i]] <- i_table$ensembl
            i_table <- subset(multidata_table, id_multidata == as.character(ligand_id[i]))
            partner_a[[i]] <- paste0("single:", i_table$name)
            receptor_a[[i]] <- i_table$receptor == 1

            if(i_table$secreted){
                flag_secreted <- TRUE
            }
        }

        # receptor
        i_table <- subset(gene_protein_table, protein_multidata_id == as.character(receptor_id[i]))
        # complex
        if(nrow(i_table) != 1){
            i_table <- subset(multidata_table, id_multidata == as.character(receptor_id[i]))
            i_receptor <- i_table$name
            partner_b[[i]] <- paste0("complex:", i_table$name)
            ensembl_b[[i]] <- ""
            receptor_b[[i]] <- i_table$receptor == 1

            if(i_table$secreted){
                flag_secreted <- TRUE
            }
            flag_integrin <- TRUE
        }else{ # single protein
            i_receptor <- i_table$gene_name
            ensembl_b[[i]] <- i_table$ensembl
            i_table <- subset(multidata_table, id_multidata == as.character(receptor_id[i]))
            partner_b[[i]] <- paste0("single:", i_table$name)
            receptor_b[[i]] <- i_table$receptor == 1

            if(i_table$secreted){
                flag_secreted <- TRUE
            }
        }

        # interaction_pair
        interaction_pair[i] <- paste0(i_ligand, "_", i_receptor)
        # secreted
        secreted[[i]] <- flag_secreted
        # is_integrin
        is_integrin[[i]] <- flag_integrin
    }

    return(data.frame(id_cp_interaction = unlist(id_cp_interaction),
                      interaction_pair = unlist(interaction_pair),
                      partner_a = unlist(partner_a),
                      partner_b = unlist(partner_b),
                      gene_a = unlist(ensembl_a),
                      gene_b = unlist(ensembl_b),
                      secreted = unlist(secreted),
                      receptor_a = unlist(receptor_a),
                      receptor_b = unlist(receptor_b),
                      annotation_strategy = filtered_interaction_table$annotation_strategy,
                      is_integrin = unlist(is_integrin),
                      check.names = FALSE))
}




#' callCellPhoneDB
#'
#' run cellphoneDB
#'
#' @param object A Seurat object
#' @param savePath
#'
#'
#' @return
#' @importFrom ggplot2 ggsave
# callCellPhoneDB <- function(object,
#                             savePath,
#                             pvalues.threshold = 0.01,
#                             means.threshold = 1.0){
#
#     sampleName <- object@project.name
#
#     savePath <- filePathCheck(savePath)
#
#     savePath <- R.utils::getAbsolutePath(savePath)
#
#     savePath_basic <- file.path(savePath, sampleName)
#     savePath_interact <- file.path(savePath_basic, "interact")
#
#     labels <- list.files(savePath_interact)
#     files <- file.path(savePath_interact, labels)
#
#     for(i in 1 : length(files)){
#         system_cmd <- "cellphonedb method statistical_analysis"
#         # meta.txt
#         system_cmd <- paste0(system_cmd, " ", file.path(files[[i]], "meta.txt"))
#         # counts.txt
#         system_cmd <- paste0(system_cmd, " ", file.path(files[[i]], "counts.txt"))
#         # project name
#         system_cmd <- paste0(system_cmd, " --project-name out")
#         # output path
#         system_cmd <- paste0(system_cmd, " --output-path ", files[[i]])
#
#         # run method
#         system(system_cmd)
#
#         # filter data
#         output.path <- file.path(files[[i]], "out")
#
#         pvalues <- read.csv(file.path(output.path, "pvalues.txt"),
#                             sep = "\t",
#                             check.names = F)
#         rownames(pvalues) <- pvalues$id_cp_interaction
#
#         means <- read.csv(file.path(output.path, "means.txt"),
#                           sep = "\t",
#                           check.names = F)
#         rownames(means) <- means$id_cp_interaction
#
#         if(!is.null(pvalues.threshold)){
#             pvalues <- pvalues[which(pvalues[[12]] <= pvalues.threshold |
#                                          pvalues[[13]] <= pvalues.threshold |
#                                          pvalues[[14]] <= pvalues.threshold |
#                                          pvalues[[15]] <= pvalues.threshold), ]
#         }
#
#         if(!is.null(means.threshold)){
#             means <- means[which(means[[12]] >= means.threshold |
#                                      means[[13]] >= means.threshold |
#                                      means[[14]] >= means.threshold |
#                                      means[[15]] >= means.threshold), ]
#         }
#
#         pvalues.rownames <- rownames(pvalues)
#         means.rownames <- rownames(means)
#
#         intersect.rownames <- intersect(pvalues.rownames, means.rownames)
#
#         pvalues <- pvalues[intersect.rownames, ]
#         means <- means[intersect.rownames, ]
#
#         write.table(pvalues,
#                     file = file.path(output.path, "filtered_pvalues.txt"),
#                     sep = "\t",
#                     quote = F,
#                     row.names = F,
#                     col.names = T)
#
#         write.table(means,
#                     file = file.path(output.path, "filtered_means.txt"),
#                     sep = "\t",
#                     quote = F,
#                     row.names = F,
#                     col.names = T)
#
#         # draw interaction area
#         meta <- read.csv(file.path(files[[i]], "meta.txt"),
#                          sep = "\t",
#                          check.names = F)
#         barcode <- meta$barcode
#         t.object <- object[, barcode]
#
#         p <- Spatial_Plot(t.object,
#                           "seurat_clusters",
#                           crop = T,
#                           discrete = T,
#                           colors = unique(t.object@meta.data[["cluster.color"]]),
#                           show.tissue = T,
#                           legend.position = "top",
#                           legend.title = "cluster",
#                           legend.size = 8)
#
#         ggsave(file.path(files[[i]], "region.png"),
#                p,
#                height = 5,
#                width = 8)
#
#         system_cmd <- "cellphonedb plot dot_plot"
#         # means.txt
#         system_cmd <- paste0(system_cmd, " --means-path ", file.path(files[[i]], "out/filtered_means.txt"))
#         # pvalues.txt
#         system_cmd <- paste0(system_cmd, " --pvalues-path ", file.path(files[[i]], "out/filtered_pvalues.txt"))
#         # output path
#         system_cmd <- paste0(system_cmd, " --output-path ", file.path(files[[i]], "out"))
#
#         # run dot plot
#         system(system_cmd)
#
#         system_cmd <- "cellphonedb plot heatmap_plot"
#         # meta.txt
#         system_cmd <- paste0(system_cmd, " ", file.path(files[[i]], "meta.txt"))
#         # pvalues.txt
#         system_cmd <- paste0(system_cmd, " --pvalues-path ", file.path(files[[i]], "out/filtered_pvalues.txt"))
#         # output path
#         system_cmd <- paste0(system_cmd, " --output-path ", file.path(files[[i]], "out"))
#
#         system(system_cmd)
#
#     }
# }
