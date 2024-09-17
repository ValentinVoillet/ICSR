#'
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Function
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Leiden algorithm
#'
#' @param data
#' @param markers Markers to use for clustering.
#' @param k Number of k neighbors.
#' @param res Resolution. A parameter controlling the coarseness of the clusters.
#' @param niter Number of iterations to run the Leiden algorithm.
#' @param seed Seed for the random number generator.
#'
#' @return partition of Leiden clusters
#' @export
#'
#' @examples
leiden_local <- function(data, markers, k = 30, res = 1, niter = 10, seed = 1234){
  #- Require
  require(igraph)
  require(RANN)
  require(Rphenograph)
  require(leiden)
  require(data.table)

  #- Data & function
  jaccard_coeff <- function(idx) {
    .Call('Rphenograph_jaccard_coeff', PACKAGE = 'Rphenograph', idx)
  }
  data_mat <- data %>%
    dplyr::select(markers) %>%
    as.matrix()

  #- Finding nearest neighbors
  snn <- RANN::nn2(data_mat, searchtype = "standard", k = k + 1)$nn.idx[, -1]

  #- Compute jaccard coefficient between nearest-neighbor sets
  links <- jaccard_coeff(snn)

  #- Build undirected graph from the weighted links
  links <- links[links[, 1] > 0, ]
  relations <- base::as.data.frame(links)
  colnames(relations) <- c("from", "to", "weight")

  #- iGraph
  g <- igraph::graph.data.frame(relations, directed = FALSE)

  #- Leiden community
  set.seed(seed)
  partition <- leiden::leiden(g, n_iterations = niter, resolution_parameter = res, seed = seed)
  return(partition)
}

