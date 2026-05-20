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
#' @param use_conda Use conda to create the python environment or not. If a path
#' is given, will use that conda env.
#' @return nothing
#' @export
#' @importFrom reticulate py_save_object
#' @importFrom Seqinfo seqinfo genome seqlevels seqlengths seqnames
#' @importFrom rtracklayer export
#' @importFrom ArchR getGenomeAnnotation getCellColData
#' @importFrom jsonlite toJSON
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
  if(length(dots$data_dir)==0){
    dots$data_dir <- outputDirectory
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
    bands <- bands[bands[, 1] %in% dots$chrom_list, ,
                   drop=FALSE]
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
  writeToJSON(dots, outputDirectory=outputDirectory)
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
    "resolution_fh" = list(500000),
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
  json <- toJSON(json, auto_unbox = TRUE, pretty = TRUE)
  writeLines(json, file.path(outputDirectory, 'Fasthigashi.JSON'))
}

#' @importFrom AnnotationHub AnnotationHub query
getCytobands <- function(genome){
  cytoband <- data.frame()
  getURL <- function(url){
    if(nchar(url)==0) return(NULL)
    headers <- curlGetHeaders(url)
    http_status <- attributes(headers)$status
    if(http_status == 200){
      response <- read.delim(text=readLines(gzcon(url(url))), header = FALSE)
      return(response)
    }else{
      return(NULL)
    }
  }
  url <- paste0('https://hgdownload.cse.ucsc.edu/goldenpath/',
                genome, '/database/cytoBand.txt.gz')
  cytoband <- getURL(url)
  if(length(cytoband)){
    return(cytoband)
  }
  ## emsembl
  ensemblURL <- function(genome){
    url <- 'https://ftp.ensembl.org/pub/'
    db <- ensembl_db[ensembl_db$assembly==genome, , drop=FALSE]
    if(nrow(db)==0){
      return('')
    }
    release <- db$release[1]
    species <- db$name[1]
    assembly <- db$assembly_ver
    return(paste0(url, 'release-', release, '/mysql/',
                  species, '_core_', release, '_',
                  assembly, '/'))
  }
  url <- ensemblURL(genome)
  # column: karyotype_id, seq_region_id, seq_region_start, seq_region_end, band, stain
  cytoband <- getURL(paste0(url, 'karyotype.txt.gz'))
  if(length(cytoband)){
    colnames(cytoband) <- c("karyotype_id", "seq_region_id",
                            "seq_region_start", "seq_region_end",
                            "band", "stain")
    # seq region column: seq_region_id, name, coord_system_id, length
    seq_region_id <- getURL(paste0(url, 'seq_region.txt.gz'))
    if(length(seq_region_id)){
      colnames(seq_region_id) <- c('seq_region_id', 'name',
                                   'coord_system_id', 'length')
      cytoband <- merge(cytoband, seq_region_id)
      cytoband <- cytoband[, c('name', 'seq_region_start', 'seq_region_end', 'band', 'stain')]
      cytoband$seq_region_start <- cytoband$seq_region_start - 1
      return(cytoband)
    }
  }
  # annotation hub
  ah = AnnotationHub()
  qy <- query(ah, c('cytoBand', genome))
  if(length(qy)>0){
    cytoband <- qy[[1]]
    cytoband <- as.data.frame(cytoband)
    cytoband$width <- NULL
    cytoband$strand <- NULL
    return(cytoband)
  }
  # biomaRt

  warning('No cytoband available in the database.')
  return(cytoband)
}

#' @importFrom httr GET content_type stop_for_status content
retreiveEnsemblInfo <- function(server = "https://rest.ensembl.org",
                                ext = "/info/species?",
                                version){

  servers <- c(
    "115"="https://sep2025.archive.ensembl.org/",
    "114"="https://may2025.archive.ensembl.org/",
    "113"="https://oct2024.archive.ensembl.org/",
    "112"="https://may2024.archive.ensembl.org/",
    "111"="https://jan2024.archive.ensembl.org/",
    "110"="https://jul2023.archive.ensembl.org/",
    "109"="https://feb2023.archive.ensembl.org/",
    "108"="https://oct2022.archive.ensembl.org/",
    "107"="https://jul2022.archive.ensembl.org/",
    "106"="https://apr2022.archive.ensembl.org/",
    "105"="https://dec2021.archive.ensembl.org/",
    "104"="https://may2021.archive.ensembl.org/",
    "103"="https://feb2021.archive.ensembl.org/",
    "102"="https://nov2020.archive.ensembl.org/",
    "101"="https://aug2020.archive.ensembl.org/",
    "100"="https://apr2020.archive.ensembl.org/",
    "99" ="https://jan2020.archive.ensembl.org/",
    "98" ="https://sep2019.archive.ensembl.org/",
    "97" ="https://jul2019.archive.ensembl.org/",
    "96" ="https://apr2019.archive.ensembl.org/",
    "95" ="https://jan2019.archive.ensembl.org/",
    "94" ="https://oct2018.archive.ensembl.org/",
    "93" ="https://jul2018.archive.ensembl.org/",
    "92" ="https://apr2018.archive.ensembl.org/",
    "91" ="https://dec2017.archive.ensembl.org/",
    "90" ="https://aug2017.archive.ensembl.org/",
    "89" ="https://may2017.archive.ensembl.org/",
    "88" ="https://mar2017.archive.ensembl.org/",
    "87" ="https://dec2016.archive.ensembl.org/",
    "86" ="https://oct2016.archive.ensembl.org/",
    "85" ="https://jul2016.archive.ensembl.org/",
    "84" ="https://mar2016.archive.ensembl.org/",
    "83" ="https://dec2015.archive.ensembl.org/",
    "82" ="https://sep2015.archive.ensembl.org/",
    "81" ="https://jul2015.archive.ensembl.org/",
    "80" ="https://may2015.archive.ensembl.org/",
    "79" ="https://mar2015.archive.ensembl.org/",
    "78" ="https://dec2014.archive.ensembl.org/",
    "77" ="https://oct2014.archive.ensembl.org/",
    "76" ="https://aug2014.archive.ensembl.org/",
    "75" ="https://feb2014.archive.ensembl.org/",
    "74" ="https://dec2013.archive.ensembl.org/",
    "73" ="https://sep2013.archive.ensembl.org/",
    "72" ="https://jun2013.archive.ensembl.org/",
    "71" ="https://apr2013.archive.ensembl.org/",
    "70" ="https://feb2013.archive.ensembl.org/",
    "69" ="https://nov2012.archive.ensembl.org/",
    "68" ="https://jul2012.archive.ensembl.org/",
    "67" ="https://may2012.archive.ensembl.org/",
    "66" ="https://feb2012.archive.ensembl.org/",
    "65" ="https://dec2011.archive.ensembl.org/",
    "64" ="https://sep2011.archive.ensembl.org/",
    "63" ="https://jun2011.archive.ensembl.org/",
    "62" ="https://apr2011.archive.ensembl.org/",
    "61" ="https://feb2011.archive.ensembl.org/",
    "60" ="https://nov2010.archive.ensembl.org/",
    "59" ="https://aug2010.archive.ensembl.org/",
    "58" ="https://may2010.archive.ensembl.org/",
    "57" ="https://mar2010.archive.ensembl.org/",
    "56" ="https://sep2009.archive.ensembl.org/",
    "55" ="https://may2009.archive.ensembl.org/",
    "54" ="https://may2009.archive.ensembl.org/"
  )

  if(!missing(version) && server=="https://rest.ensembl.org"){
    version <- as.character(version)
    if(version %in% names(servers)){
      server <- sub('archive', 'rest', servers[version])
    }else{
      stop('server name not available!')
    }
  }
  r <- GET(paste(server, ext, sep = ""), content_type("application/json"))

  stop_for_status(r)

  content <- content(r)[['species']]

  n <- unique(unlist(lapply(content, names)))

  db <- do.call(rbind, lapply(content, function(.ele){
    vapply(.ele[n], function(x) paste(x, collapse=';'), character(1L))
  }))

  db <- data.frame(db)
}

