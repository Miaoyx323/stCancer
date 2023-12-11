SpatialStructureDetection <- function(
        object,
        feature,
        slice = 'slice1',
        global.Moran.I = 0.5,
        p.adjust.method = 'fdr',
        p.thres.local = 0.05,
        p.adj.thres.local = 0.05,
        p.thres.z.score = 0.05,
        median.filter = TRUE,
        z.score.filter = TRUE
){
    require(Seurat)
    require(SeuratObject)
    
    pos <- GetTissueCoordinates(object[[slice]])
    pos.dist <- dist(x = pos)
    pos.dist.mat <- as.matrix(x = pos.dist)
    
    weights <- 1 / pos.dist.mat^2
    
    diag(x = weights) <- 0
    pos.dist.mat.min <- min(pos.dist.mat[pos.dist.mat > 0])
    weights[pos.dist.mat > pos.dist.mat.min * 2.5] <- 0
    
    ROWSUM <- rowSums(weights)
    ROWSUM[ROWSUM == 0] <- 1
    weights <- weights / ROWSUM
    
    X <- scale(object[[feature]])[, 1]
    
    global.res <- RunGlobalMoransI(X = X,
                                   weights = weights)
    if(global.res$I < global.Moran.I){
        print(paste0('I: ', global.res$I,
                     ', less than the global threshold ', global.Moran.I))
        return(object)
    }
    
    print(paste0('I: ', global.res$I))
    
    local.res <- RunLocalMoransI(X = X,
                                 weights = weights,
                                 p.adjust.method = p.adjust.method)
    local.res['area.MoransI'] <- ((local.res$P.local < p.thres.local) &
                                      (local.res$I > 0) & (local.res$ZI > 0))
    if(median.filter){
        local.res['area.MoransI'] <- local.res['area.MoransI'] & (X > median(X))
    }
    
    if(z.score.filter){
        local.res['area.filter'] <- FALSE
        tmp.object <- object[, local.res['area.MoransI'] == TRUE]
        tmp.X <- scale(tmp.object[[feature]])[, 1]
        local.res[names(tmp.X), 'area.filter'] <- RunZscorePvalue(
            tmp.X,
            alternative = 'greater'
        ) < p.thres.z.score
    }else{
        local.res['area.filter'] <- FALSE
        local.res['area.filter'] <- (local.res$area.MoransI & 
                                         (local.res$P.adj < p.adj.thres.local))
    }
    
    colnames(local.res) <- paste0(feature, '.', colnames(local.res))
    object <- AddMetaData(object,
                          local.res)
    return(object)
}

RunGlobalMoransI <- function(
        X,
        weights
){
    require(ape)
    res <- ape::Moran.I(X, weights)
    res$ZI <- (res$observed - res$expected) / res$sd
    res <- as.data.frame(res)
    res <- res[c("observed","expected","sd","ZI","p.value")]
    colnames(res) <- c('I', 'EI', 'SdI', 'ZI', 'P')
    return(res)
}

RunLocalMoransI <- function(
        X,
        weights,
        p.adjust.method = 'fdr'
){
    require(spdep)
    res <- data.frame(matrix(nrow = length(X),
                             ncol = 5))
    rownames(res) <- names(X)
    colnames(res) <- c('I', 'EI', 'VarI', 'ZI', 'P.local')
    
    # remove isolated areas
    ROWSUM <- rowSums(weights)
    X.use <- X[ROWSUM != 0]
    weights.use <- weights[ROWSUM != 0, ROWSUM != 0]
    
    suppressMessages(suppressWarnings(
        listw <- mat2listw(weights.use)
    ))
    
    tmp <- localmoran(X.use, listw)
    colnames(tmp) <- c('I', 'EI', 'VarI', 'ZI', 'P.local')
    
    res[rownames(tmp), ] <- tmp
    res[which(rowSums(is.na(res)) > 0), ] <- c(0,0,0,0,1)
    
    res <- as.data.frame(res)
    res$P.adj <- p.adjust(res$P,
                          method = p.adjust.method)
    return(res)
}

RunZscorePvalue <- function(
        X,
        alternative = 'two.sided'
){
    X <- scale(X,
               center = TRUE,
               scale = TRUE)
    if (alternative == "two.sided") {
        pv <- 2 * pnorm(abs(X), lower.tail = FALSE)
    }
    else if (alternative == "greater") {
        pv <- pnorm(X, lower.tail = FALSE)
    }
    else {
        pv <- pnorm(X)
    }
    return(pv)
}