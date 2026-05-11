#' import RNA reads and DNA reads as a ArchR project
#' @param rnaFiles,atacFiles A named characters.
#' The file path of RNA and DNA count matrix files.
#' @param bmap,barcodesFun The parameter required by \link{replace_RNA_barcode}.
#' @param outputDirectory The output folder the the merged RNA and DNA
#' ArchR project
#' @param annotationGenes If there is not standard annotations available,
#' please input the gene annotations in GRanges object with symbols metadata
#' @param verbose Print message or not.
#' @return An ArchRProject object.
#' @export
#' @importFrom ArchR getArchRGenome getGeneAnnotation import10xFeatureMatrix createArrowFiles ArchRProject subsetArchRProject saveArchRProject getCellNames addGeneExpressionMatrix
#' @importFrom Seurat CreateSeuratObject Read10X_h5 Read10X
#' @importFrom GenomicFeatures genes
#' @importFrom SummarizedExperiment rowRanges<-
createArchRproj <- function(
    rnaFiles, atacFiles,
    bmap, barcodesFun=function(x){
      do.call(rbind, strsplit(colnames(x), '_'))[, 4]
    },
    outputDirectory = "ArchR_proj",
    annotationGenes = NULL,
    verbose=TRUE
    ){
  stopifnot(length(rnaFiles)==length(atacFiles))
  stopifnot(length(rnaFiles)>0)
  stopifnot(all(names(atacFiles) %in% names(rnaFiles)))
  stopifnot(all(names(rnaFiles) %in% names(atacFiles)))
  rnaFiles <- rnaFiles[names(atacFiles)]
  barLen <- 18
  checkBmap(bmap, barLen = barLen)
  if(is.null(annotationGenes)){
    currentGenome <- getArchRGenome()
    if(is.null(currentGenome)){
      stop('The genome is not available. Please use addArchRGenome to set genome.')
    }
    geneAnno <- getGeneAnnotation()
    annotationGenes <- geneAnno$genes
  }

  stopifnot(is(annotationGenes, 'GRanges'))
  stopifnot('symbol is not available in metadata of annotationGenes'=
              length(annotationGenes)==length(annotationGenes$symbol))

  ## RNA
  if(verbose) message('Reading RNA reads')
  data <- lapply(rnaFiles, function(.ele){
    if(grepl('.h5$', .ele)){
      Read10X_h5(.ele)
    }else{
      Read10X(dirname(.ele))
    }
  })

  rna <- lapply(seq_along(rnaFiles), function(i){
    seu <- CreateSeuratObject(counts = data[[i]], project = names(rnaFiles)[i], min.cells=0, min.features=0)#donot filter out any cells !important
    seu$sample <- names(rnaFiles)[i]
    seu
  })

  if(verbose) message('replace the RNA barcodes')
  rna <- replace_RNA_barcode(rna, bmap=bmap, barcodesFun=barcodesFun)

  if(verbose) message('write the RNA reads to a tmp file')
  tmpF <- tempfile(pattern = names(rnaFiles), fileext = "h5")
  null <- mapply(Write10X_h5, rna, tmpF)

  if(verbose) message('import the tmp file by import10xFeatureMatrix')
  seRNA <- import10xFeatureMatrix(
    input = tmpF,
    names = names(rnaFiles),
    strictMatch = TRUE,
    force = TRUE
  )

  if(verbose) message('fix the rowRanges information in seRNA')
  rn <- rownames(seRNA)
  seRNA <- seRNA[rn %in% genes$symbol]
  rn <- rownames(seRNA)
  rr <- genes[match(rn, genes$symbol)]
  names(rr) <- rn
  rowRanges(seRNA) <- rr

  if(verbose) message('create arrow files')
  ArrowFiles <- createArrowFiles(
    inputFiles = atacFiles,
    sampleNames = names(atacFiles),
    outputNames = file.path(tempdir(), names(atacFiles)),
    minTSS = 1,
    minFrags = 1000,
    addTileMat = TRUE,
    addGeneScoreMat = TRUE
  )

  od <- file.path(tempdir(), 'ATAC_ArchR')
  if(verbose) message('create ATAC ArchR project in ', od)
  dir.create(od, recursive = TRUE, showWarnings = FALSE)
  dna <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = od, copyArrows = TRUE)

  if(verbose){
    message('There are ', length(getCellNames(dna)), ' cells in DNA reads. ',
            'There are ', length(colnames(seRNA)), ' cells in RNA reads. ',
            'The shared barcodes number is ',
            length(which(getCellNames(dna) %in% colnames(seRNA))), ' (',
            round(100*length(which(getCellNames(dna) %in%
                                     colnames(seRNA)))/
                    length(getCellNames(dna)), digits = 2),
            '%).' )
  }

  if(verbose) message('Writing ArchR files to ', outputDirectory)
  cellsToKeep <- which(getCellNames(dna) %in% colnames(seRNA))
  proj <- subsetArchRProject(ArchRProj = dna, cells = getCellNames(dna)[cellsToKeep], outputDirectory = outputDirectory, force = TRUE)

  if(verbose) message('Add gene expression matrix to the project')
  proj <- addGeneExpressionMatrix(input = proj, seRNA = seRNA, strictMatch = TRUE, force = TRUE)

  proj <- saveArchRProject(ArchRProj = proj, overwrite = TRUE, load = TRUE)

  return(proj)
}
