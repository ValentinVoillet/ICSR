#' Leiden Clustering with Optional CHOIR Pruning
#'
#' Constructs a Shared Nearest Neighbor (SNN) graph from pre-normalized marker
#' data, runs the Leiden community detection algorithm, and prunes over-clustered
#' partitions using CHOIR's random forest statistical test.
#'
#' @param data       A data frame of cells × features (pre-normalized).
#' @param markers    Character vector of column names to use as features.
#' @param k          Number of nearest neighbors for SNN graph (default: 30).
#' @param res        Leiden resolution parameter — higher = more clusters (default: 2).
#' @param niter      Number of Leiden iterations (default: 10).
#' @param seed       Random seed for reproducibility (default: 1234).
#' @param n_cores    Number of cores for CHOIR's pruneTree (default: 4).
#'
#' @return A named list:
#'   \item{clusters}{Raw Leiden cluster labels (character vector).}
#'   \item{pruned_clusters}{CHOIR-pruned cluster labels (character vector).}
#'   \item{graph}{igraph network.}
#'
#' @import igraph RANN leiden
#' @importFrom Rphenograph Rphenograph
#' @importFrom dplyr select all_of n_distinct
#' @importFrom Matrix sparseMatrix
#' @importFrom SingleCellExperiment SingleCellExperiment reducedDim
#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors metadata
#' @importFrom CHOIR inferTree pruneTree
#'
#' @examples
#' \dontrun{
#'   clusters <- leiden_CHOIR_local(data = filtered_exprs, markers = my_markers, k = 30, res = 0.8)
#' }
#'
#' @export
leiden_CHOIR_local <- function(data,
                               markers,
                               k       = 30,
                               res     = 2,
                               niter   = 10,
                               seed    = 1234,
                               n_cores = 4,
                               min.accuracy = 0.5) {

  # ── 1. Feature matrix ───────────────────────────────────────────────────────
  # Extract only the marker columns and convert to a plain numeric matrix.
  # All downstream steps (RANN, igraph, CHOIR) require a matrix, not a data frame.
  mat <- data %>%
    dplyr::select(dplyr::all_of(markers)) %>%
    as.matrix()

  cell_ids <- paste0("cell_", seq_len(nrow(mat)))
  rownames(mat) <- cell_ids


  # ── 2. k-Nearest Neighbor search ────────────────────────────────────────────
  # RANN::nn2 finds the k+1 nearest neighbors for each cell (Euclidean distance).
  # We request k+1 because the closest "neighbor" is always the cell itself
  # (distance = 0), so we drop the first column to get k true neighbors.
  message("Finding nearest neighbors (k = ", k, ")...")
  nn_idx <- RANN::nn2(mat, searchtype = "standard", k = k + 1)$nn.idx[, -1]


  # ── 3. Jaccard similarity + SNN graph ───────────────────────────────────────
  # For each pair of cells that share neighbors, compute the Jaccard coefficient:
  #   J(A, B) = |neighbors(A) ∩ neighbors(B)| / |neighbors(A) ∪ neighbors(B)|
  # This rewards pairs that share many of the *same* neighbors, not just proximity.
  # The Rphenograph C++ backend makes this fast for large cell counts.
  message("Computing Jaccard coefficients...")
  jaccard_coeff <- function(idx) {
    .Call('Rphenograph_jaccard_coeff', PACKAGE = 'Rphenograph', idx)
  }

  links <- jaccard_coeff(nn_idx)
  links <- links[links[, 1] > 0, ]   # drop zero-weight edges (no shared neighbors)

  relations <- as.data.frame(links)
  colnames(relations) <- c("from", "to", "weight")

  # Build the undirected weighted SNN graph
  g <- igraph::graph_from_data_frame(relations, directed = FALSE)


  # ── 4. Leiden community detection ───────────────────────────────────────────
  # Leiden improves on Louvain by guaranteeing well-connected communities.
  # Higher `res` → more, smaller clusters. With res = 2 this intentionally
  # over-clusters, relying on CHOIR below to merge indistinct clusters.
  message("Running Leiden algorithm (res = ", res, ")...")
  set.seed(seed)
  partition <- leiden::leiden(g,
                              resolution_parameter = res,
                              n_iterations         = niter,
                              seed                 = seed)
  clusters  <- as.character(partition)
  message(sprintf("  → %d raw Leiden clusters", dplyr::n_distinct(clusters)))


  # ── 5. Binary NN and SNN sparse matrices ────────────────────────────────────
  # CHOIR's pruneTree() needs:
  #   nn_matrix  — binary cell × cell matrix (1 = cell j is a k-NN of cell i)
  #   snn_matrix — weighted SNN adjacency from the Leiden graph
  n_cells    <- nrow(mat)
  nn_matrix  <- Matrix::sparseMatrix(
    i    = rep(seq_len(n_cells), each = k),
    j    = as.vector(nn_idx),
    x    = 1,
    dims = c(n_cells, n_cells)
  )
  rownames(nn_matrix) <- cell_ids
  colnames(nn_matrix) <- cell_ids
  snn_matrix <- igraph::as_adjacency_matrix(g, attr = "weight", sparse = TRUE)
  rownames(snn_matrix) <- cell_ids
  colnames(snn_matrix) <- cell_ids


  # ── 6. Build SingleCellExperiment for CHOIR ─────────────────────────────────
  # Data is already normalized upstream — stored as 'normcounts'.
  # normalization_method = "none" in pruneTree() ensures CHOIR does not attempt
  # to re-normalize, which would corrupt already-normalized values.
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(normcounts = t(mat))   # SCE convention: features × cells
  )
  colnames(sce) <- cell_ids

  # Store the marker matrix as a 'PCA'-slot reduction.
  # inferTree() and pruneTree() use this to compute centroid distances between clusters.
  SingleCellExperiment::reducedDim(sce, "PCA") <- mat   # cells × features
  rownames(SingleCellExperiment::reducedDim(sce, "PCA")) <- cell_ids


  # ── 7. inferTree(): build a hierarchical tree over Leiden clusters ───────────
  # CHOIR requires a binary merge tree (dendrogram) over flat cluster labels.
  # inferTree() computes cluster centroids in the embedding space, then builds
  # an agglomerative hierarchy based on centroid distances — clusters that are
  # most similar are merged last (root) vs. first (leaves near the tips).
  # This tree structure is what pruneTree() traverses bottom-up.
  message("Building CHOIR hierarchy from Leiden clusters...")
  named_clusters <- setNames(clusters, cell_ids)

  cluster_tree <- CHOIR::inferTree(
    cluster_labels = named_clusters,
    reduction      = mat              # cells × features
  )


  # ── 8. Inject tree into SCE metadata ────────────────────────────────────────
  # CHOIR reads everything from metadata(sce)$CHOIR.
  # We must pre-populate this slot so pruneTree() can find the tree.
  S4Vectors::metadata(sce) <- list(
    CHOIR = list(
      clusters = list(cluster_tree = cluster_tree)
    )
  )


  # ── 9. pruneTree(): statistically merge indistinct clusters ─────────────────
  # pruneTree() walks the tree bottom-up. At each node it trains a random forest
  # to classify cells from the two child clusters, then runs a permutation test:
  #   H0: the two clusters are indistinguishable (RF accuracy ≈ chance)
  # If the test fails to reject H0 (p > alpha), the two clusters are merged.
  # This continues up the tree until all remaining splits are statistically supported.
  #
  # Key parameters:
  #   alpha        — FDR threshold; 0.05 is standard
  #   use_assay    — 'normcounts' matches the pre-normalized input stored in step 6
  #   nn_matrix    — used to identify local neighborhoods for RF train/test split
  #   snn_matrix   — the Leiden SNN graph weights
  #   reduction    — embedding used for centroid distances (consistent with inferTree)
  message("Pruning tree with CHOIR random forest tests...")
  set.seed(seed)
  sce_pruned <- CHOIR::pruneTree(
    object               = sce,
    key                  = "CHOIR",
    alpha                = 0.05,
    input_matrix         = t(mat),
    use_assay            = "normcounts",   # pre-normalized; no re-normalization needed
    cluster_tree         = cluster_tree,
    nn_matrix            = nn_matrix,
    snn_matrix           = snn_matrix,
    reduction            = mat,
    normalization_method = "none",
    random_seed          = seed,
    verbose              = TRUE,
    n_cores              = n_cores,
    min_accuracy         = min.accuracy
  )


  # ── 10. Extract final cluster labels ────────────────────────────────────────
  pruned_clusters <- as.character(
    SummarizedExperiment::colData(sce_pruned)$CHOIR_clusters_0.05
  )

  message(sprintf(
    "Leiden: %d clusters  →  CHOIR pruned: %d clusters",
    dplyr::n_distinct(clusters),
    dplyr::n_distinct(pruned_clusters)
  ))


  # ── 11. Return ───────────────────────────────────────────────────────────────
  return(list(
    clusters        = clusters,          # raw Leiden labels
    pruned_clusters = pruned_clusters,   # CHOIR-validated labels
    graph = g
  ))
}
