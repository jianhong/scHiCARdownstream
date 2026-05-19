#' create contact matrix npy file for scDeepLUCIA
#' @param proj An ArchRProject object.
#' @param pairs_files File paths for pairs.
#' The input pairs file must have columnes barcode, chr1, pos1, chr2, pos2
#' @param outdir The output folder for contact npy files.
#' @param groupBy The column names of cell column data used for split the cells
#' @param chroms The chromosome should be worked on.
#' @param resolution The resolution for the genomic feature npy file.
#' @param ... The parameters not used.
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
#' @importFrom dplyr count %>% .data
#' @importFrom GenomicRanges tileGenome
#' @importFrom BSgenome getSeq
#' @examples
#' # example code
#'
createContactMatrix <- function(
    proj,
    pairs_files,
    outdir='scDeepLUCIA_inputs/sliced_con_array',
    groupBy='Sample',
    chroms, resolution=5000,
    ...,
    verbose = TRUE, use_conda = FALSE){
  np <- import('numpy')
  stopifnot(is(proj, 'ArchRProject'))
  stopifnot(length(names(pairs_files))==length(pairs_files))
  anno <- getGenomeAnnotation(proj)
  prefix <- 'contact_030M.'
  genome <- get(anno$genome)
  seqlev <- standardChromosomes(genome)
  if(!missing(chroms)){
    chroms <- chroms[chroms %in% seqlev]
  }else{
    chroms <- seqlev
  }
  stopifnot('No chromosome detectable.'=length(chroms)>0)

  outdir <- file.path(outdir, genome(genome)[chroms[1]])
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  ccd <- getCellColData(proj)
  ccd <- ccd[, colnames(ccd) %in% groupBy, drop=FALSE]
  if(ncol(ccd)!=length(groupBy)){
    stop('Not all groupBy available in CellColData of input project.')
  }

  barcodes <- getArchRBarcode(proj)
  barcodes <- split(barcodes, make.names(apply(ccd, 1, paste, sep='_')))
  if(!all(names(pairs) %in% names(barcodes))){
    bgp <- paste(names(barcodes), collapse=', ')
    stop('The names of pairs are not all in the barcodes groups. ',
         'The names of barcodes groups are: ', bgp)
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
  })
  if(verbose) message('Filter the barcodes')
  pairs <- lapply(pairs, function(.ele){
    .ele[.ele[, 1] %in% barcodes, ]
  })
  pairs <- split(pairs, names(pairs))

  pairs <- mapply(function(pair, bc){
    pair <- do.call(rbind, pair)
    rownames(pair) <- NULL
    pair[pair[, 1] %in% bc, ]
  }, pairs, barcodes[names(pairs)], SIMPLIFY = FALSE)

  if(verbose) message('Write the unsorted pair to tsv file')
  null <- mapply(function(ps, n){
    pff <- file.path(outdir, n)
    dir.create(pff, showWarnings = FALSE)
    write.table(ps, file = file.path(pff, paste0(n, '.contact.pairs')),
                quote = FALSE, sep = '\t',
                col.names = FALSE, row.names = FALSE)
    if(verbose){
      file_name <- file.path(pff, paste0(n, '.contact.pairs'))
      message('command lines: ',
              'bgzip ', file_name)
      srt_file <- paste0(sub('.contact', '.sort', file_name), '.gz')
      message('pairtools sort -o ', srt_file,
              ' ', file_name, '.gz')
      message('pairix -f ', srt_file)
      cool_file <- sub('.contact.pairs', '.', resolution,'.cool', file_name)
      message('cooler cload pairix --max-split 2 --nproc ${ncore} ${chrom_size_file}:',
              resolution, ' ',
              srt_file, ' ', cool_file)
      message('cooler zoomify --balance -r ', resolution, 'N -n ${ncore} -o ',
              sub('.cool', '.mcool', cool_file), ' ', cool_file)
    }
  },pairs, names(pairs), SIMPLIFY = FALSE)

  seqlen <- seqlengths(genome)[chroms]
  gr <- tileGenome(seqlen, tilewidth = resolution,
                   cut.last.tile.in.chrom = TRUE)
  gr <- split(gr, seqnames(gr))
  l <- lengths(gr)

  null <- mapply(function(ps, n){
    pff <- file.path(outdir, n)
    colnames(ps) <- c('barcodes', 'seq1', 'start1', 'seq2', 'start2')
    ps[, 'start1'] <- ceiling(ps[, 'start1']/resolution)
    ps[, 'start2'] <- ceiling(ps[, 'start2']/resolution)
    ps <- ps[, -1, drop=FALSE]
    ps <- ps[ps$seq1==ps$seq2, , drop=FALSE]
    counts_df <- ps %>% count(.data$seq1, .data$start1, .data$start2)
    counts_df <- split(counts_df[, c('start1', 'start2', 'n')], counts_df[, 'seq1'])
    for(chr in names(counts_df)){
      mx <- matrix(data=0, nrow=l[chr], ncol=l[chr])
      for(i in seq.int(nrow(counts_df[[chr]]))){
        mx[counts_df[[chr]][i, 'start1'], counts_df[[chr]][i, 'start2']] <-
          counts_df[[chr]][i, 'n']
      }
      np$save(file.path(pff, paste0(prefix, chr, ".npy")), mx)
      rm(mx)
    }
  },pairs, names(pairs), SIMPLIFY = FALSE)
}
