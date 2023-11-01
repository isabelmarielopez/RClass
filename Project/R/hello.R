#' readSNV
#'
#' This function takes one argument, which is the raw SNV data. The function returns a SEURAT object from the data, which is simply the matrix processed into a format friendly to R.
#'
#' @param path The path to the raw SNV data
#'
#' @return The Seurat object representing the data
#'
#' @export

readSNV <- function(matrixData) {
  # ...
}

readSNV <- function(matrixData){
  colnames(matrixData) = gsub('-1', paste0('_', "a"), colnames(matrixData))
  a <- CreateSeuratObject(matrixData, min.cells =5)
  return(a)
}

#' Seur.QC
#'
#' This function takes one argument, which is the seurat object representing the SNV data. The function returns a plot visualizing the genes per nucleus of the sample.
#'
#' @param obj The seurat object created by readSNV
#'
#' @return modified object
#'
#' @export

Seur.QC <- function(object) {
  # ...
}

Seur.QC <- function(obj){
  object <- Add_Mito_Ribo_Seurat(seurat_object = obj, species = "mouse")
  object <- Add_Cell_Complexity_Seurat(seurat_object = object)
  object$mitoRatio = object@meta.data$percent_mito/100
  object$log10GenesPerUMI = log10(object$nFeature_RNA)/log10(object$nCount_RNA)
  object$nUMI = object$nCount_RNA
  object$log10_nUMI = log10(object$nUMI)
  object$condition = object$orig.ident
  object$tissue = "Nonlesion" #variable
  object$treatment = "PBS" #variable
  object$nGene = object$nFeature_RNA
  object@assays$RNA@meta.features$cpn = Matrix::rowMeans(object@assays$RNA@counts)
  return(object)
}

#' plots
#'
#' This function takes one argument, which is the seurat object representing the SNV data. The function returns a UMAP and gene by nuclei plot of the sample.
#'
#' @param object The seurat object created by readSNV
#'
#' @export
plots <- function(object) {
  # ...
}

plots <- function(object){
  object<-Seur.QC(object)
  object = SCTransform(object, vst.flavor = "v2", variable.features.n = 3000)
  object = RunPCA(object, verbose = TRUE, assay = "SCT", npcs = 100)
  object = FindNeighbors(object, dims = 1:40, verbose = FALSE, reduction = "pca")
  object = FindClusters(object, resolution = 1.5, verbose = FALSE) #variables = resolution
  object = RunUMAP(object, dims = 1:40, n.neighbors = 30L, min.dist = 0.45, verbose = TRUE, assay = "SCT") #variables = n.neibors, dims, and min.dist
  DimPlot_scCustom(object, figure_plot = TRUE, label.size = 11, repel = TRUE, pt.size = 1.0, reduction = "umap")
}
