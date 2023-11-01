#' readSNV
#'
#' This function takes one argument, which is the raw SNV data. The function returns a SEURAT object from the data, which is simply the matrix processed into a format friendly to R.
#'
#' @name readSNV
#' 
#' @param matrix The path to the raw SNV data
#'
#' @return The Seurat object representing the data
#'
#' @export


readSNV <- function(matrixData){
  colnames(matrixData) = gsub('-1', paste0('_', "a"), colnames(matrixData))
  a <- CreateSeuratObject(matrixData, min.cells =5)
  return(a)
}


