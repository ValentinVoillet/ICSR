#' Run UMAP Dimensionality Reduction
#'
#' This function is a wrapper for the \code{uwot::umap} implementation. It
#' transforms high-dimensional flow cytometry data into a low-dimensional
#' space (typically 2D) for visualization.
#'
#' @param dt A \code{data.table} or \code{data.frame} containing the single-cell events.
#' @param markers Character vector. The names of the markers (columns) to use for the projection.
#' @param n_neighbors Integer. The size of local neighborhood used for manifold approximation.
#'   Larger values result in more global structure. Default is 10.
#' @param min_dist Numeric. The effective minimum distance between embedded points.
#'   Smaller values result in tighter clusters. Default is 0.1.
#' @param n_components Integer. The dimension of the space to embed into. Default is 2.
#' @param verbose Logical. If \code{TRUE}, logs progress details to the console.
#' @param seed Integer. Seed for reproducibility. Default is 1234.
#'
#' @return A \code{matrix} of UMAP embeddings with \code{n_components} columns.
#'
#' @import uwot data.table dplyr
#' @export
#'
#' @examples
#' # umap_coords <- UMAP_local(dt = dt.exprs, markers = clustering_markers, n_neighbors = 15)
UMAP_local <- function(dt, markers, n_neighbors = 10, min_dist = .1, verbose = TRUE, n_components = 2, seed = 1234){
  # 1. Setup and Matrix Conversion
  # ---------------------------------------------------------------------------
  # UMAP requires a numeric matrix input
  data_mat <- dt %>%
    dplyr::select(dplyr::all_of(markers)) %>%
    as.matrix()

  # 2. Execute UMAP Algorithm
  # ---------------------------------------------------------------------------
  if(verbose) message("Running UMAP dimensionality reduction...")
  set.seed(seed)
  umap_embedding <- uwot::umap(
    X = data_mat,
    n_neighbors = n_neighbors,
    min_dist = min_dist,
    n_components = n_components,
    verbose = verbose
  )

  # 3. Output
  # ---------------------------------------------------------------------------
  return(umap_embedding)
}

