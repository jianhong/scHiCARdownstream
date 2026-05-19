#' create One-hot encoded genomic feature npy file.
#' @param proj The ArchR project
#' @param outdir The output directory
#' @param chroms The chromosome should be worked on.
#' @param resolution The resolution for the genomic feature npy file.
#' @return nothing
#' @importFrom GenomeInfoDb standardChromosomes
#' @importFrom Seqinfo seqlengths genome
#' @importFrom reticulate import
#' @export
#' @examples
#' # example code
#'
createOnehotGenomicNpy <- function(
    proj, outdir='scDeepLUCIA_inputs/sliced_epi_array',
    chroms, resolution=5000){
  stopifnot(is(proj, 'ArchRProject'))
  anno <- getGenomeAnnotation(proj)
  prefix <- 'isHC-seq_onehot.'
  genome <- get(anno$genome)
  seqlev <- standardChromosomes(genome)
  if(!missing(chroms)){
    chroms <- chroms[chroms %in% seqlev]
  }else{
    chroms <- seqlev
  }
  stopifnot('No chromosome detectable.'=length(chroms)>0)
  outdir <- file.path(outdir, genome(genome)[chroms[1]], 'isHC')
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  seqlen <- seqlengths(genome)[chroms]
  gr <- tileGenome(seqlen, tilewidth = resolution,
                   cut.last.tile.in.chrom = TRUE)
  gr <- split(gr, seqnames(gr))

  MAPPING = rbind(
    "A"= c(1, 0, 0, 0),
    "C"= c(0, 1, 0, 0),
    "G"= c(0, 0, 1, 0),
    "T"= c(0, 0, 0, 1),
    "N"= c(1, 1, 1, 1),
    "a"= c(1, 0, 0, 0),
    "c"= c(0, 1, 0, 0),
    "g"= c(0, 0, 1, 0),
    "t"= c(0, 0, 0, 1),
    "n"= c(1, 1, 1, 1)
  )
  mode(MAPPING) <- 'integer'
  np <- import('numpy')
  for(n in names(gr)){ ## use for loop to handle the chromosome one by one
    seq <- getSeq(genome, gr[[n]])
    encode <- lapply(seq, function(reads){
      MAPPING[strsplit(as.character(reads), '')[[1]], ]
    })
    encode[[length(encode)]] <-
      rbind(encode[[length(encode)]],
            matrix(0L, nrow=resolution, ncol=4))[seq.int(resolution), ]
    encode <- simplify2array(encode)
    encode <- aperm(encode, c(3, 1, 2))
    np$save(file.path(outdir, paste0(prefix, n, ".npy")), encode)
    rm(encode, seq)
  }
}
