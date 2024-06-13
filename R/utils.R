#' getDefaultColors
#'
#' Get default colors if not required
#'
#' @return A list of colors
#'
#'
getDefaultColors <- function(n = NULL, type = 1){
    if(type == 1){
        colors <- c("#ff1a1a", "#1aff1a", "#1a1aff", "#ffff1a", "#ff1aff",
                    "#ff8d1a", "#7cd5c8", "#c49a3f", "#5d8d9c", "#90353b",
                    "#507d41", "#502e71", "#1B9E77", "#c5383c", "#0081d1",
                    "#674c2a", "#c8b693", "#aed688", "#f6a97a", "#c6a5cc",
                    "#798234", "#6b42c8", "#cf4c8b", "#666666", "#ffd900",
                    "#feb308", "#cb7c77", "#68d359", "#6a7dc9", "#c9d73d")
    }else if(type == 2){
        if(n <= 14){
            colors <- c("#437BFE", "#FEC643", "#43FE69", "#FE6943", "#C643FE",
                        "#43D9FE", "#B87A3D", "#679966", "#993333", "#7F6699",
                        "#E78AC3", "#333399", "#A6D854", "#E5C494")
        }
        else if(n <= 20){
            colors <- c("#87b3d4", "#d5492f", "#6bd155", "#683ec2", "#c9d754",
                        "#d04dc7", "#81d8ae", "#d34a76", "#607d3a", "#6d76cb",
                        "#ce9d3f", "#81357a", "#d3c3a4", "#3c2f5a", "#b96f49",
                        "#4e857e", "#6e282c", "#d293c8", "#393a2a", "#997579")
        }else if(n <= 30){
            colors <- c("#628bac", "#ceda3f", "#7e39c9", "#72d852", "#d849cc",
                        "#5e8f37", "#5956c8", "#cfa53f", "#392766", "#c7da8b",
                        "#8d378c", "#68d9a3", "#dd3e34", "#8ed4d5", "#d84787",
                        "#498770", "#c581d3", "#d27333", "#6680cb", "#83662e",
                        "#cab7da", "#364627", "#d16263", "#2d384d", "#e0b495",
                        "#4b272a", "#919071", "#7b3860", "#843028", "#bb7d91")
        }else{
            colors <- c("#982f29", "#5ddb53", "#8b35d6", "#a9e047", "#4836be",
                        "#e0dc33", "#d248d5", "#61a338", "#9765e5", "#69df96",
                        "#7f3095", "#d0d56a", "#371c6b", "#cfa738", "#5066d1",
                        "#e08930", "#6a8bd3", "#da4f1e", "#83e6d6", "#df4341",
                        "#6ebad4", "#e34c75", "#50975f", "#d548a4", "#badb97",
                        "#b377cf", "#899140", "#564d8b", "#ddb67f", "#292344",
                        "#d0cdb8", "#421b28", "#5eae99", "#a03259", "#406024",
                        "#e598d7", "#343b20", "#bbb5d9", "#975223", "#576e8b",
                        "#d97f5e", "#253e44", "#de959b", "#417265", "#712b5b",
                        "#8c6d30", "#a56c95", "#5f3121", "#8f846e", "#8f5b5c")
        }
    }else if(type == 3){
        colors <- c("#588dd5", "#c05050", "#07a2a4", "#f5994e",
                    "#9a7fd1", "#59678c", "#c9ab00", "#7eb00a")
    }else if(type == 4){
        colors <- c("#FC8D62", "#66C2A5", "#8DA0CB", "#E78AC3",
                    "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")
    }else if(type == 5){
        colors <- c("#c14089", "#6f5553", "#E5C494", "#738f4c",
                    "#bb6240", "#66C2A5", "#2dfd29", "#0c0fdc")
    }
    if(!is.null(n)){
        if(n <= length(colors)){
            colors <- colors[1:n]
        }else{
            step <- 16777200 %/% (n - length(colors)) - 2
            add.colors <- paste0("#", as.hexmode(seq(from = sample(1:step, 1),
                                                     by = step, length.out = (n-length(colors)))))
            colors <- c(colors, add.colors)
        }
    }
    return(colors)
}

#' getDefaultDimColors
#'
#' Return the colors for discrete features
#'
#' @import RColorBrewer
#'
getDefaultDimColors <- function(n){
    if(n <= 11){
        # brewer set3
        colors <- c('#8DD3C7','#FFFFB3','#BEBADA','#FB8072','#80B1D3','#FDB462',
                    '#B3DE69','#FCCDE5','#BC80BD','#CCEBC5','#FFED6F')
    }else if(n <= 25){
        # poly wright25
        colors <- c('#1F78C8','#FF0000','#33A02C','#6A33C2','#FF7F00','#565656',
                    '#FFD700','#A6CEE3','#FB6496','#B2DF8A','#CAB2D6','#FDBF6F',
                    '#999999','#EEE685','#C8308C','#FF83FA','#C814FA','#0000FF',
                    '#36648B','#00E2E5','#00FF00','#778B00','#BEBE00','#8B3B00',
                    '#A52A3C')
    }else{
        colors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(n)
    }
    return(colors[1:n])
}

#' getDefaultFeatureColors
#'
#' Return the colors for continuous features
#'
#' @import RColorBrewer
#'
getDefaultFeatureColors <- function(type = 'seq',
                                    name = NULL){
    if(type == 'seq'){
        colors <- colorRampPalette(
            colors = brewer.pal(n = 11, name = ifelse(is.null(name),
                                                      'YlOrRd',
                                                      name))
        )
    }else if(type == 'div'){
        colors <- colorRampPalette(
            colors = rev(x = brewer.pal(n = 11, name = ifelse(is.null(name),
                                                              'RdYlBu',
                                                              name)))
        )
    }
    return(colors)
}


#' getGeneManifest
#'
#' Read genes.tsv
#'
#' @return A dataframe of gene information
#'
#' @importFrom utils read.delim
#'
#'
getGeneManifest <- function(data.path, h5=FALSE, filter.matrix=TRUE){
    if(!h5){
        data.path <- file.path(data.path, "filtered_feature_bc_matrix")
        gene.loc <- file.path(data.path, 'features.tsv.gz')
        if(!grepl("feature_bc_matrix", data.path)) {
            gene.loc <- file.path(data.path, 'genes.tsv')
        }
        gene.manifest <- read.delim(gene.loc, header = FALSE, stringsAsFactors = FALSE)
        colnames(gene.manifest) <-
            c("EnsemblID", "Symbol", "Type")[1:dim(gene.manifest)[2]]
    }else{
        matrix <- Read10X_h5(file.path(data.path, ifelse(filter.matrix,
                                                         'filtered_feature_bc_matrix.h5',
                                                         'raw_feature_bc_matrix.h5')))
        gene.manifest <- data.frame(Symbol = rownames(matrix))
    }

    return(gene.manifest)
}

#' addGeneAnno
#'
#' Annotate mito & ribo genes
#'
#' @return Gene manifest
#'
addGeneAnno <- function (gene.manifest,
                         species = "human"){

    annotation <- rep("other", dim(gene.manifest)[1])

    if(species == "human"){
        mito.genes <- grep('^MT-', rownames(gene.manifest), value = TRUE)
        ribo.genes <- grep('^RPL|^RPS|^MRPL|^MRPS', rownames(gene.manifest), value = TRUE)
        diss.genes <- read.table(file.path(system.file(package = 'stCancer'), 'txt/diss-genes.txt'),
                                 header = F, stringsAsFactors = F)$V1
    }else if(species == "mouse"){
        mito.genes <- grep('^mt-', rownames(gene.manifest), value = TRUE)
        ribo.genes <- grep('^Rpl|^Rps|^Mrpl|^Mrps', rownames(gene.manifest), value = TRUE)
        diss.genes <- read.table(file.path(system.file(package = 'stCancer'), 'txt/diss-genes.txt'),
                                 header = F, stringsAsFactors = F)$V1
        diss.genes <- getMouseGene(diss.genes)
    }

    annotation[rownames(gene.manifest) %in% mito.genes] <- "mitochondrial"
    annotation[rownames(gene.manifest) %in% ribo.genes] <- "ribosome"
    annotation[rownames(gene.manifest) %in% diss.genes] <- "dissociation"

    gene.manifest$Annotation <- annotation
    return(gene.manifest)
}


#' limitData
#'
#' limit data between min and max
#'
#' @param data data
#' @param min min
#' @param max max
#'
#' @return A new data
#'
#'
limitData <- function(data,
                      min = NULL,
                      max = NULL){
    data2 <- data
    if(!is.null(min)){
        data2[data2 < min] <- min
    }
    if(!is.null(max)){
        data2[data2 > max] <- max
    }
    return(data2)
}


#' combine_spot
#'
#' Combine spot.manifest and other matrix
#'
#' @param spot.manifest
#'
#' @return Spot manifest
#'
combine_spot <- function(spot.manifest,
                         mat){
    rownames(spot.manifest) <- spot.manifest$barcode

    mat.colname <- colnames(mat)

    for(i in 1 : length(mat.colname)){
        spot.manifest[[mat.colname[i]]] <- rep(0, nrow(spot.manifest))
        spot.manifest[rownames(mat), mat.colname[i]] <- mat[, i]
    }

    return(spot.manifest)
}


#' #' filePathCheck
#' #'
#' #' check the format of the file path
#' #'
#' #' @param filePath
#' #'
#' #' @return filePath
#' #'
#' #'
#' filePathCheck <- function(filePath){
#'     return(gsub("\\\\", "/", filePath))
#' }

#' getMouseGene
#'
#' @import homologene
#'
#'
getMouseGene <- function(hg.genes){
    res <- homologene(hg.genes, inTax = 9606, outTax = 10090)
    return(res[['10090']])
}

# getMouseGene <- function(hg.genes, bool.name = F, deduplicate = T){
#     hg.mm.HomologyGenes <- read.table(system.file("txt", "hg-mm-HomologyGenes.txt", package = "stCancer"),
#                                       header = T, stringsAsFactors = F)
#     hg.mm.HomologyGenes <- subset(hg.mm.HomologyGenes, hgGenes %in% hg.genes)
#
#     if(deduplicate){
#         hg.num <- table(hg.mm.HomologyGenes$hgGenes)
#         hg.mm.HomologyGenes <- subset(hg.mm.HomologyGenes, !(hgGenes %in% names(hg.num)[hg.num > 1]))
#         mm.num <- table(hg.mm.HomologyGenes$mmGenes)
#         hg.mm.HomologyGenes <- subset(hg.mm.HomologyGenes, !(mmGenes %in% names(mm.num)[mm.num > 1]))
#     }
#
#     mm.genes <- hg.mm.HomologyGenes$mmGenes
#
#     if(bool.name){
#         names(mm.genes) <- hg.mm.HomologyGenes$hgGenes
#     }
#     return(mm.genes)
# }

get_neighbors_indices <- function(row, col){
    return(list(c(row-1, col-1),
                c(row-1, col+1),
                c(row, col-2),
                c(row, col+2),
                c(row+1, col-1),
                c(row+1, col+1)))
}
