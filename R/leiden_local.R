#' Leiden-based Clustering for Single-Cell Data
#'
#' This function performs unsupervised clustering using the Leiden algorithm
#' on a Shared Nearest Neighbor (SNN) graph. It follows the Rphenograph
#' logic by calculating Jaccard coefficients but uses the more robust Leiden
#' community detection to find cell populations.
#'
#' @param data A \code{data.table} or \code{data.frame} containing single-cell events.
#' @param markers Character vector. The names of the markers (columns) to use for clustering.
#' @param k Integer. The number of nearest neighbors used to build the SNN graph. Default is 30.
#' @param res Numeric. The resolution parameter controlling the granularity of the
#'   clusters. Higher values lead to more clusters. Default is 1.
#' @param niter Integer. The number of iterations to run the Leiden algorithm. Default is 10.
#' @param seed Integer. Seed for reproducibility. Default is 1234.
#'
#' @return A character vector of cluster assignments (partitions) for each cell.
#'
#' @details
#' The function utilizes the \code{Rphenograph} Jaccard coefficient calculation
#' (via C++ call) to weigh the edges of the graph before passing it to the
#' Leiden algorithm.
#'
#' @import igraph RANN Rphenograph leiden data.table
#' @export
#'
#' @examples
#' # clusters <- leiden_local(data = filtered_exprs, markers = my_markers, k = 30, res = 0.8)
leiden_local <- function(data, markers, k = 30, res = 1, niter = 10, seed = 1234){
  # 1. Setup and Matrix Conversion
  # ---------------------------------------------------------------------------
  # The Jaccard coefficient calculation is optimized for matrices
  data_mat <- data %>%
    dplyr::select(dplyr::all_of(markers)) %>%
    as.matrix()

  # 2. Shared Nearest Neighbor (SNN) Graph Construction
  # ---------------------------------------------------------------------------
  message("Finding nearest neighbors (k = ", k, ")...")
  # Find k+1 neighbors because the first neighbor is always the point itself
  snn <- RANN::nn2(data_mat, searchtype = "standard", k = k + 1)$nn.idx[, -1]
  message("Computing Jaccard coefficients...")
  # Use the Rphenograph C++ backend for fast Jaccard coefficient calculation
  jaccard_coeff <- function(idx) {
    .Call('Rphenograph_jaccard_coeff', PACKAGE = 'Rphenograph', idx)
  }
  links <- jaccard_coeff(snn)
  # Filter for meaningful links (weight > 0)
  links <- links[links[, 1] > 0, ]
  relations <- base::as.data.frame(links)
  colnames(relations) <- c("from", "to", "weight")

  # 3. Graph and Community Detection
  # ---------------------------------------------------------------------------
  # Build an undirected graph where weights represent shared neighbor density
  g <- igraph::graph_from_data_frame(relations, directed = FALSE)
  message("Running Leiden algorithm (res = ", res, ")...")
  set.seed(seed)
  partition <- leiden::leiden(
    g,
    n_iterations = niter,
    resolution_parameter = res,
    seed = seed
  )
  # Return the cluster assignments as a factor/character
  return(as.character(partition))
}

