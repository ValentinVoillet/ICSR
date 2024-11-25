#'
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Function
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' UMAP
#'
#' @param dt
#' @param markers Markers to use for UMAP.
#' @param n_neighbors The size of local neighborhood. See `uwot::umap`.
#' @param min_dist The effective minimum distance between embedded points. See `uwot::umap`.
#' @param n_components The dimension of the space to embed into. See `uwot::umap`.
#' @param verbose If TRUE, log details to the console. See `uwot::umap`.
#' @param seed Seed for the random number generator.
#'
#' @return UMAP embeddings
#' @export
#'
#' @examples
UMAP_local <- function(dt, markers, n_neighbors = 10, min_dist = .1, verbose = TRUE, n_components = 2, seed = 1234){
  #- Require
  require(uwot)
  require(data.table)

  #- Data
  data_mat <- dt %>%
    dplyr::select(all_of(markers)) %>%
    as.matrix()

  #- UMAP
  set.seed(seed)
  umap <- uwot::umap(data_mat, n_neighbors = n_neighbors, min_dist = min_dist, n_components = n_components, verbose = verbose)

  #- Output
  return(umap)
}

