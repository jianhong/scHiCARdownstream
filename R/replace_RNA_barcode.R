#' Replace RNA barcodes to matched ATAC barcodes
#' @param rna A list of Seurat Object.
#' @param bmap A named character vector. The names of the vector are
#' RNA barcode and the elements are DNA barcode.
#' @param barcodesFun How to extract the RNA barcode.
#' The default function is used to extract RNA barcodes for the original manuscript. The barcodes is the cell names in the format of "xw329_day0_120_TACAGCAAACTCAACGTA". If your cell names is in the format of "TACAGCAAACTCAACGTA-1", you may want to try:
#' barcodesFun = function(x)\{do.call(rbind, strsplit(colnames(x), '-)\[1\])\}
#' @return A list of Seurat object with metadata DNAbarcode as replaced barcode
#' @importFrom Biostrings reverseComplement DNAStringSet
#' @export
#' @examples
#' # example code
#'
replace_RNA_barcode <- function(
    rna, bmap,
    barcodesFun=function(x){
      do.call(rbind, strsplit(colnames(x), '_'))[, 4]
    }){
  barLen <- 18
  checkBmap(bmap, barLen = barLen)
  stopifnot(is.function(barcodesFun))
  # split barcode to 1,2,3
  rna <- lapply(rna, function(.ele){
    stopifnot('rna must be a list of Seurat object'=
                is(.ele, 'Seurat'))
    ## remove the
    barcodes <- barcodesFun(.ele)
    stopifnot('all barcodes must be a character of length 18'=
                all(nchar(barcodes)==barLen))
    barcodes <- do.call(rbind, strsplit(barcodes, paste0("(?<=.{", barLen/3, "})"), perl=TRUE))
    for(i in seq.int(3)){
      .ele[[paste0('barcode', i)]] <- barcodes[, i]
    }
    .ele
  })

  rna <- lapply(rna, function(.ele){
    tobereplaced <- sum(.ele$barcode1 %in% names(bmap))
    b1 <- .ele$barcode1
    if(tobereplaced>0){
      k <- b1 %in% names(bmap)
      b1[k] <- bmap[b1[k]]
    }
    b2 <- .ele$barcode2
    b3 <- as.character(reverseComplement(DNAStringSet(.ele$barcode3)))
    .ele$DNAbarcode <- paste0(b1, b2, b3)
    .ele
  })

  return(rna)
}

checkBmap <- function(bmap, barLen=18){
  barLen <- barLen/3
  stopifnot(is.character(bmap))
  stopifnot(length(names(bmap))==length(bmap))
  stopifnot(all(nchar(bmap)==barLen))
  stopifnot(all(nchar(names(bmap))==barLen))
}
