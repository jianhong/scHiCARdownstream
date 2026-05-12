#' create files required for FastHigashi
#' @param proj An ArchRProject object.
#' @param pairs_files File paths for pairs.
#' The input pairs file must have columnes barcode, chr1, pos1, chr2, pos2
#' @param outputDirectory The output folder for filtered pairs and config files
#' of Higashi.
#' @param labelInfoColumns The column names of cell column data used for
#' label_info.pickle
#' @param plot_label,... The parameters for config JSON.
#' See https://github.com/ma-compbio/Higashi/wiki/Configure-the-parameters-for-Fast-Higashi.
#' @param verbose Print message or not.
#' @param use_conda Use conda to create the python enviroment or not
#' @return nothing
#' @export
#' @importFrom reticulate py_save_object
#' @importFrom karyoploteR getCytobands
#' @importFrom Seqinfo seqinfo genome seqlevels seqlengths seqnames
#' @importFrom rtracklayer export
#' @importFrom ArchR getGenomeAnnotation getCellColData
#' @importFrom rjson toJSON
#' @importFrom utils write.table read.delim
#' @examples
#' # example code
#'
prepareHigashiFiles <- function(
    proj,
    pairs_files,
    outputDirectory='per_cell',
    labelInfoColumns=c('Sample', 'Clusters_Combined'),
    plot_label='Clusters_Combined',
    ...,
    verbose = TRUE, use_conda = FALSE){
  stopifnot(is(proj, 'ArchRProject'))
  if(!dir.exists(outputDirectory)) {
    dir.create(outputDirectory, recursive = TRUE)
  }
  ccd <- getCellColData(proj)
  ccd <- ccd[, colnames(ccd) %in% labelInfoColumns, drop=FALSE]
  if(ncol(ccd)!=length(labelInfoColumns)){
    stop('Not all labelInfoColumns available in CellColData of input project.')
  }
  barcodes <- getArchRBarcode(proj)
  barLen <- nchar(barcodes[1])
  genes <- getGeneAnnotation(proj)[['genes']]
  chrom_list <- seqlevels(genes)
  chrom_list <- chrom_list[!grepl("Y", chrom_list, ignore.case = TRUE)]
  dots <- list(...)
  if(length(dots$chrom_list)==0){
    dots$chrom_list <- chrom_list
  }
  if(length(dots$impute_list)==0){
    dots$impute_list <- chrom_list
  }
  if(length(dots$genome_reference_path)==0){
    chrom.size <- seqlengths(genes)[chrom_list]
    chrom.size <- as.data.frame(chrom.size)
    write.table(chrom.size, file.path(outputDirectory, 'chrom.sizes.txt'),
                quote = FALSE, col.names = FALSE, sep = '\t')
    dots$genome_reference_path <- file.path(outputDirectory, 'chrom.sizes.txt')
  }
  if(length(dots$cytoband_path)==0){
    genome <- genome(seqinfo(genes))[1]
    bands <- getCytobands(genome=genome)
    bands <- bands[as.character(seqnames(bands)) %in% dots$chrom_list]
    bands <- as.data.frame(bands)
    bands$width <- NULL
    bands$strand <- NULL
    write.table(bands, file.path(outputDirectory, 'cytoBand.txt'),
                quote = FALSE, col.names=FALSE, sep = '\t')
    dots$cytoband_path <- file.path(outputDirectory, 'cytoBand.txt')
  }
  if(length(dots$batch_id)==0){
    dots$batch_id <- 'batch'
    if(!dots$batch_id %in% colnames(ccd)){
      ccd[[dots$batch_id]] <- 'batch_1'
    }
  }
  if(length(dots$library_id)==0){
    dots$library_id <- 'Sample'
  }
  if(length(dots$blacklist)==0){
    anno <- getGenomeAnnotation(proj)
    blacklist <- anno$blacklist
    if(length(blacklist)){
      export(blacklist, file.path(outputDirectory, 'blacklist.bed'))
      dots$blacklist <- file.path(outputDirectory, 'blacklist.bed')
    }
  }

  if(verbose) message('Read Pairs')
  pairs <- lapply(pairs_files, read.delim, header=FALSE)
  null <- lapply(pairs, function(.ele){
    stopifnot('The input pairs file must have columnes barcode, chr1, pos1, chr2, pos2'=
                is.integer(.ele[, 3]))
    stopifnot('The input pairs file must have columnes barcode, chr1, pos1, chr2, pos2'=
                is.integer(.ele[, 5]))
    stopifnot('The input pairs file must have columnes barcode, chr1, pos1, chr2, pos2'=
                is.character(.ele[, 1]))
    if(nchar(.ele[1, 1])!=barLen){
      stop('The barcode in input pairs file does not have ', barLen,
           ' characters.')
    }
  })
  if(verbose) message('Filter the barcodes')
  pairs <- lapply(pairs, function(.ele){
    .ele <- .ele[.ele[, 1] %in% barcodes, ]
    split(.ele[, -1], .ele[, 1])
  })
  if(verbose) message('Save the data')
  filelist <- mapply(function(psl, n){
    pff <- file.path(outputDirectory, n)
    dir.create(pff, showWarnings = FALSE)
    mapply(function(ps, cellname){
      psn <- file.path(pff, paste0(cellname, '.tsv'))
      write.table(ps, file = psn, quote = FALSE, sep = '\t',
                  col.names = FALSE, row.names = FALSE)
      return(psn)
    }, psl, names(psl))
  },pairs, names(pairs), SIMPLIFY = FALSE)
  barcodes2 <- unlist(lapply(pairs, names), use.names = FALSE)
  filelist <- unlist(filelist)
  label_info <- ccd[match(barcodes2, barcodes), , drop=FALSE]
  writeLines(filelist, file.path(outputDirectory, 'filelist.txt'))
  writeToPickle(label_info, outputDirectory = outputDirectory,
                use_conda=use_conda)
  writeToJSON(, outputDirectory=outputDirectory)
}

getArchRBarcode <- function(proj){
  barcodes <- gsub(".*#", "", getCellNames(proj))
}

writeToPickle <- function(label_info, outputDirectory, use_conda=FALSE){
  if(use_conda){
    prepare_miniconda()
    prepare_env()
  }
  label_info_list <- as.list(label_info)
  py_save_object(label_info_list,
                 file.path(outputDirectory, "label_info.pickle"))
}

writeToJSON <- function(additional=list(), outputDirectory){
  json <- list(
    "config_name" = "scHiCAR",
    "data_dir" = ".",
    "temp_dir" = tempdir(),
    "genome_reference_path" = "./chrom.sizes.txt",
    "cytoband_path" = "./cytoBand.txt",
    "chrom_list" = additional$chrom_list,
    "resolution" = 1000000,
    "resolution_cell" = 1000000,
    "resolution_fh" = 500000,
    "minimum_distance" = 2000000,
    "maximum_distance" = -1,
    "local_transfer_range" = 1,
    "dimensions" = 64,
    "impute_list" = additional$impute_list,
    "minimum_impute_distance" = 0,
    "maximum_impute_distance" = -1,
    "neighbor_num" = 5,
    "plot_start" = 0,
    "plot_end" = -1,
    "plot_label" = "Clusters_Combined",
    "batch_id" = "batch",
    "call_tads" = FALSE,
    "embedding_name" = "exp_zinb3",
    "embedding_epoch" = 80,
    "no_nbr_epoch" = 80,
    "with_nbr_epoch" = 60,
    "cpu_num" = -1,
    "gpu_num" = 1,
    "optional_smooth" = FALSE,
    "optional_quantile" = FALSE,
    "rank_thres" = 1,
    "loss_mode" = "zinb",
    "random_walk" = FALSE,
    "UMAP_params" = list("n_neighbors" = 30, "min_dist" = 0.3),
    "TSNE_params" = list("n_neighbors" = 15),
    "input_format" ="higashi_v2",
    "header_included" = FALSE,
    "contact_header" = c("chrom1", "pos1", "chrom2", "pos2"),
    "structured" = TRUE
  )
  for(i in names(additional)){
    json[[i]] <- additional[[i]]
  }
  json <- toJSON(json)
  writeLines(json, file.path(outputDirectory, 'Fasthigashi.JSON'))
}
