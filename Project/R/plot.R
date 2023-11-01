#' plots
#'
#' This function takes one argument, which is the seurat object representing the SNV data. The function returns a UMAP and gene by nuclei plot of the sample.
#' @name plots
#'
#' @param object The seurat object created by readSNV
#'
#' @export

plots <- function(object){
  object<-Seur.QC(object)
  object = SCTransform(object, vst.flavor = "v2", variable.features.n = 3000)
  object = RunPCA(object, verbose = TRUE, assay = "SCT", npcs = 100)
  object = FindNeighbors(object, dims = 1:40, verbose = FALSE, reduction = "pca")
  object = FindClusters(object, resolution = 1.5, verbose = FALSE) #variables = resolution
  object = RunUMAP(object, dims = 1:40, n.neighbors = 30L, min.dist = 0.45, verbose = TRUE, assay = "SCT") #variables = n.neibors, dims, and min.dist
  DimPlot_scCustom(object, figure_plot = TRUE, label.size = 11, repel = TRUE, pt.size = 1.0, reduction = "umap")
}