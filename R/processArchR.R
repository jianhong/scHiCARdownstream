#' process the ATAC/RNA merged ArchR project by standard pipeline
#' @param proj An ArchRProject object. The output of \link{createArchRproj}.
#' @param resolutions The resolutions for each steps.
#' @param minDist A number that determines how tightly the UMAP is allowed to pack points together. This argument is passed to min_dist in uwot::umap(). For more info on this see https://jlmelville.github.io/uwot/abparams.html.
#' @param skip_rmDup skip filter doublets or not.
#' @param verbose Print message or not.
#' @return An ArchRProject object.
#' @export
#' @importFrom ArchR addDoubletScores filterDoublets addIterativeLSI addCombinedDims addUMAP addClusters saveArchRProject
#' @examples
#' # example code
#'
processArchR <- function(
    proj,
    resolutions=c("LSI_ATAC"=0.2, "LSI_RNA"=0.2,
                  "Clusters_ATAC"=0.4, "Clusters_RNA"=0.4,
                  'Clusters_Combined'=0.4),
    minDist = 0.8,
    skip_rmDup = FALSE,
    verbose=TRUE){
  stopifnot(is(proj, 'ArchRProject'))
  stopifnot(all(c('LSI_ATAC', 'LSI_RNA', 'Clusters_ATAC', 'Clusters_RNA', 'Clusters_Combined') %in% names(resolutions)))
  stopifnot(is.numeric(resolutions))
  stopifnot(is.numeric(minDist))

  if(!skip_rmDup){
    if(verbose) message('Add doublet score to the project and filter doublets')
    proj <- addDoubletScores(proj, force = TRUE)
    proj <- filterDoublets(proj)
  }

  if(verbose) message('Add LSI matrix for ATAC')
  proj <- addIterativeLSI(
    ArchRProj = proj,
    clusterParams = list(
      resolution = resolutions[["LSI_ATAC"]],
      sampleCells = 10000,
      n.start = 10
    ),
    saveIterations = FALSE,
    useMatrix = "TileMatrix",
    depthCol = "nFrags",
    name = "LSI_ATAC"
  )
  if(verbose) message('Add LSI matrix for RNA')
  proj <- addIterativeLSI(
    ArchRProj = proj,
    clusterParams = list(
      resolution = resolutions[["LSI_RNA"]],
      sampleCells = 10000,
      n.start = 10
    ),
    saveIterations = FALSE,
    useMatrix = "GeneExpressionMatrix",
    depthCol = "Gex_nUMI",
    varFeatures = 2500,
    firstSelection = "variable",
    binarize = FALSE,
    name = "LSI_RNA"
  )
  proj <- addCombinedDims(proj,
                          reducedDims = c("LSI_ATAC", "LSI_RNA"),
                          name =  "LSI_Combined")
  proj <- addUMAP(proj,
                  reducedDims = "LSI_ATAC",
                  name = "UMAP_ATAC",
                  minDist = minDist,
                  force = TRUE)
  proj <- addUMAP(proj,
                  reducedDims = "LSI_RNA",
                  name = "UMAP_RNA",
                  minDist = minDist,
                  force = TRUE)
  proj <- addUMAP(proj,
                  reducedDims = "LSI_Combined",
                  name = "UMAP_Combined",
                  minDist = minDist,
                  force = TRUE)
  proj <- addClusters(proj,
                      reducedDims = "LSI_ATAC",
                      name = "Clusters_ATAC",
                      resolution = resolutions[['Clusters_ATAC']],
                      force = TRUE)
  proj <- addClusters(proj,
                      reducedDims = "LSI_RNA",
                      name = "Clusters_RNA",
                      resolution = resolutions[['Clusters_RNA']],
                      force = TRUE)
  proj <- addClusters(proj,
                      reducedDims = "LSI_Combined",
                      name = "Clusters_Combined",
                      resolution = resolutions[['Clusters_Combined']],
                      force = TRUE)

  proj <- saveArchRProject(ArchRProj = proj, overwrite = TRUE, load = TRUE)

  proj
}
