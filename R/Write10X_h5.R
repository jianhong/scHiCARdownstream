#' Write 10X hdf5 file
#'
#' Write count matrix to 10X CellRanger hdf5 file.
#' This file will be used for `ArchR::import10xFeatureMatrix` function.
#'
#' @param seu The Seurat object after barcode remapping by \link{replace_RNA_barcode}.
#' @param filename Path to h5 file
#' @param layer,assay The parameter for \link[SeuratObject:AssayData]{GetAssayData}.
#' @param barcodes The colname of barcode in metadata. If it is a Function, such as colnames, The function will be used to extract the  barcode. Otherwise, the column in metadata will be used as barcode.
#' @importFrom hdf5r H5File
#' @importFrom Seurat GetAssayData
#' @return nothing
#' @export
#' @importFrom methods is
#' @examples
#' # example code
#'
#'
Write10X_h5 <- function(seu, filename, layer='counts', assay='RNA', barcodes='DNAbarcode') {
  if (!requireNamespace('hdf5r', quietly = TRUE)) {
    stop("Please install hdf5r to read HDF5 files")
  }
  # message("filename ", filename)
  outfile <- hdf5r::H5File$new(filename = filename, mode = 'w')
  on.exit(outfile$close_all())
  # create groupfeature_slot
  genome <- outfile$create_group("matrix")
  data <- GetAssayData(seu, layer=layer, assay=assay)
  genome[["data"]] <- as.integer(data@x)
  genome[["indices"]] <- data@i
  if (is(data, 'dgCMatrix64')) {
    genome[["indptr"]] <- as.numeric(data@p)
  } else {
    genome[["indptr"]] <- as.integer(data@p)
  }
  genome[["shape"]] <- as.integer(data@Dim)
  if(is.function(barcodes)){
    genome[["barcodes"]] <- barcodes(seu)
  }else{
    genome[["barcodes"]] <- seu[[barcodes]]
  }
  features <- genome$create_group('features')
  features[['name']] <- rownames(data)
  features[['feature_type']] <- rep('Gene Expression', nrow(data))
  outfile$close_all()
  on.exit()
}
