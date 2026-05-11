#' call peaks
#' @param proj An ArchRProject object. The output of \link{processArchR}.
#' @param groupBy The groups used to call peaks by function \link[ArchR:addReproduciblePeakSet]{addReproduciblePeakSet}.
#' @param pathToMacs2 The macs2 path.
#' @param ... Other parameters used by \link[ArchR:addReproduciblePeakSet]{addReproduciblePeakSet}.
#' @param reducedDims The reductions used to link the gene and peaks.
#' @param verbose Print message or not.
#' @return An ArchRProject object.
#' @importFrom ArchR addGroupCoverages addReproduciblePeakSet addPeakMatrix addPeak2GeneLinks findMacs2
#' @export
#' @examples
#' # example code
#'
callPeaks <- function(
    proj,
    groupBy='Clusters_Combined',
    pathToMacs2=findMacs2(),
    ...,
    reducedDims='LSI_Combined',
    verbose=TRUE){
  stopifnot(is(proj, 'ArchRProject'))
  proj <- addGroupCoverages(ArchRProj = proj, groupBy = groupBy, verbose = FALSE)
  proj <- addReproduciblePeakSet(ArchRProj = proj, groupBy = groupBy, pathToMacs2 = pathToMacs2, ...)
  proj <- addPeakMatrix(ArchRProj = proj)
  ## could we use hic links?
  proj <- addPeak2GeneLinks(ArchRProj = proj,
                            reducedDims = reducedDims,
                            useMatrix = "GeneExpressionMatrix")

  proj <- saveArchRProject(ArchRProj = proj, overwrite = TRUE, load = TRUE)

  proj
}
