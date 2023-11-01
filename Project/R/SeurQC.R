
#' Seur.QC
#'
#' This function takes one argument, which is the seurat object representing the SNV data. The function returns a plot visualizing the genes per nucleus of the sample.
#' @name Seur.QC
#'
#' @param obj The seurat object created by readSNV
#'
#' @return modified object
#'
#' @export


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
