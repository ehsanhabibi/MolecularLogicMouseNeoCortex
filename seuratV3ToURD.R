#' Import URD from Seurat
#' 
#' If you're previously performed analyses in Seurat, you can copy them over
#' directly into an URD object for analyzing pseudotime, pseudolineage, and
#' building a specification tree.
#' 
#' @param seurat.object A Seurat object
#' @return An URD object
#' @export

seuratV3ToURD <- function(seurat.object) {
  if (requireNamespace("Seurat", quietly = TRUE)) {
    # Create an empty URD object
    ds <- new("URD")
    
    # Copy over data
    #mat <- GetAssayData(object = seurat.object)
    #ds@logupx.data <- as(as.matrix(mat), "dgCMatrix")
    ds@logupx.data <- GetAssayData(object = seurat.object)
    if(!any(dim(GetAssayData(object = seurat.object, slot = "counts")) == 0)) ds@count.data <- GetAssayData(object = seurat.object, slot = "counts")[rownames(GetAssayData(object = seurat.object)), colnames(GetAssayData(object = seurat.object))]
    
    # Copy over metadata
    ## TO DO - grab kmeans clustering info
    get.data <- NULL
    if (.hasSlot(seurat.object, "data.info")) { 
      get.data <- as.data.frame(seurat.object@data.info)
    } else if (.hasSlot(seurat.object, "meta.data")) { 
      get.data <- as.data.frame(seurat.object@meta.data) 
    }
    if(!is.null(get.data)) {
      di <- colnames(get.data)
      m <- grep("res|cluster|Res|Cluster", di, value=T, invert = T) # Put as metadata if it's not the result of a clustering.
      discrete <- apply(get.data, 2, function(x) length(unique(x)) / length(x))
      gi <- di[which(discrete <= 0.015)]
      ds@meta <- get.data[,m,drop=F]
      ds@group.ids <- get.data[,gi,drop=F]
    }
    
    # Copy over var.genes
    if(length(VariableFeatures(object = seurat.object) > 0)) ds@var.genes <- VariableFeatures(object = seurat.object)
    
    # Move over tSNE projection
    if (.hasSlot(seurat.object, "tsne.rot")) {
      if(!any(dim(seurat.object@tsne.rot) == 0)) {
        ds@tsne.y <- as.data.frame(seurat.object@tsne.rot)
        colnames(ds@tsne.y) <- c("tSNE1", "tSNE2")
      }
    #} else if (.hasSlot(seurat.object, "dr")) {
    } else{
      if(("tsne" %in% names(seurat.object)) && !any(dim(seurat.object[["tsne"]]) == 0)) {
        ds@tsne.y <- as.data.frame(seurat.object[["tsne"]]@cell.embeddings)
        colnames(ds@tsne.y) <- c("tSNE1", "tSNE2")
      }
    }
    
    # Move over PCA results
    if (.hasSlot(seurat.object, "pca.x")) {
      if(!any(dim(seurat.object@pca.x) == 0)) {
        ds@pca.load <- seurat.object@pca.x
        ds@pca.scores <- seurat.object@pca.rot
        warning("Need to set which PCs are significant in @pca.sig")
      }
      ## TO DO: Convert SVD to sdev
    #} else if (.hasSlot(seurat.object, "dr")) {
    } else{ 
      if(("pca" %in% names(seurat.object)) && !any(dim(seurat.object[["pca"]]@feature.loadings) == 0)) {
        ds@pca.load <- as.data.frame(seurat.object[["pca"]]@feature.loadings)
        ds@pca.scores <- as.data.frame(seurat.object[["pca"]]@cell.embeddings)
        ds@pca.sdev <- seurat.object[["pca"]]@stdev  #seurat.object@dr$pca@sdev
        ds@pca.sig <- pcaMarchenkoPastur(M=dim(ds@pca.scores)[1], N=dim(ds@pca.load)[1], pca.sdev=ds@pca.sdev)
      }
    }
    return(ds)
  } else {
    stop("Package Seurat is required for this function. To install: install.packages('Seurat')\n")
  }
}
