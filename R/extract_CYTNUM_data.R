#' Extract Cytokine Polyfunctionality Counts
#'
#' This function calculates the number of cytokines expressed per cell (\code{CYTNUM})
#' and summarizes the total number of cytokine-positive vs. cytokine-negative cells
#' within a defined parent population.
#'
#' @param gs A \code{GatingSet} object.
#' @param parent_node Character. The parent gate to filter cells (e.g., "CD4"). Only cells in this gate are returned.
#' @param cytokine_nodes Character vector. Subset of \code{output_nodes} used to count positive cytokines per cell.
#' @param pData_cols Character vector. Columns from \code{pData(gs)} to include in the output.
#'   Default includes BATCH, PTID, STIM, VISITNO, Run Num, Replicate, SAMP_ORD
#' @param stim_to_exclude Character vector. List of stimulations to filter out from the final results.
#'
#' @return A \code{data.table} summarized by sample and positivity status,
#'   containing total cell counts (\code{NSUB}) and cytokine-positive counts (\code{CYTNUM}).
#'
#' @import flowWorkspace tidyverse data.table
#' @export
extract_CYTNUM_data <- function(gs,
                                parent_node,
                                cytokine_nodes,
                                pData_cols = c("BATCH", "PTID", "STIM", "VISITNO", "Run Num", "Collection Num", "Replicate", "SAMP_ORD"),
                                stim_to_exclude = NULL)
{
  # 1. Main Extraction Loop (Per Sample)
  # ---------------------------------------------------------------------------
  exprs.tmp <- flowWorkspace::lapply(gs, function(x)
  {
    # --- Boolean Positivity Call ---
    # Retrieve indices for the parent and all cytokine gates
    marker_response <- NULL
    target_nodes <- c(parent_node, cytokine_nodes)
    # Retry loop to handle potential I/O or memory errors during index retrieval
    while (is.null(marker_response)) {
      marker_response <- tryCatch({
        res <- lapply(target_nodes, function(mrkr) flowWorkspace::gh_pop_get_indices(x, mrkr))
        names(res) <- target_nodes
        dplyr::bind_rows(res)
      }, error = function(e) NULL)
    }

    # --- Merging Metadata & Results ---
    # Extract requested pData columns
    pd_subset <- flowWorkspace::pData(x) %>%
      dplyr::select(dplyr::any_of(pData_cols)) %>%
      dplyr::rename_with(~ "RUNNUM", dplyr::matches("Run Num|Collection Num")) %>%
      dplyr::rename_with(~ "REPLICATE", dplyr::matches("Replicate"))
    pd_subset$FCS <- rownames(pd_subset)
    standardized_pd_names <- colnames(pd_subset)
    # Merging
    dt.res <- marker_response %>%
      dplyr::bind_cols(pd_subset[rep(1, flowWorkspace::gh_pop_get_stats(x, "root")$count), ])

    # --- Cytokine Calculation ---
    # Calculate how many cytokines each individual cell expresses
    dt.res$CYTNUM <- rowSums(dt.res[, cytokine_nodes, with = FALSE] == TRUE)

    # --- Population Summarization ---
    # Filter for parent cells and calculate totals
    dt.output <- dt.res %>%
      dplyr::filter(get(parent_node) == TRUE) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(standardized_pd_names))) %>%
      dplyr::mutate(NSUB = n()) %>%
      dplyr::mutate(boolean_CYTNUM = factor(x = ifelse(CYTNUM >= 1, TRUE, FALSE), levels = c(TRUE, FALSE))) %>%
      dplyr::group_by(NSUB, boolean_CYTNUM, .add = TRUE, .drop = FALSE) %>%
      dplyr::summarize(CYTNUM = n(), .groups = "drop")
    return(dt.output)
  })

  # 2. Final Filtering and Aggregation
  # ---------------------------------------------------------------------------
  res_final <- dplyr::bind_rows(exprs.tmp)
  if (!is.null(stim_to_exclude)) {
    stim_col <- intersect(c("STIM", "Stim"), colnames(res_final))[1]
    if (!is.na(stim_col)) {
      res_final <- res_final %>% dplyr::filter(!(get(stim_col) %in% stim_to_exclude))
    }
  }
  return(data.table::as.data.table(res_final))
}

