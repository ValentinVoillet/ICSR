#' Extract Flow Cytometry Expression Data and Population Membership
#'
#' This function iterates through a GatingSet to extract single-cell expression data
#' (Fluorescence Intensities) and boolean gating results. It supports multiple
#' transformation types (Biexponential, Arcsinh, and Asymmetric Arcsinh).
#'
#' @param gs A \code{GatingSet} object.
#' @param output_nodes Character vector. Names of the FlowJo gates to extract boolean positivity (e.g., "IFNg+").
#' @param parent_node Character. The parent gate to filter cells (e.g., "CD4"). Only cells in this gate are returned.
#' @param cytokine_nodes Character vector. Subset of \code{output_nodes} used to count positive cytokines per cell.
#' @param pData_cols Character vector. Columns from \code{pData(gs)} to include in the output.
#'   Default includes BATCH, PTID, STIM, VISITNO, Run Num, Replicate, SAMP_ORD
#' @param do.comp Logical. If \code{TRUE}, extracts compensated raw values.
#' @param do.biexp Logical. If \code{TRUE}, extracts FlowJo's biexponential transformed values.
#' @param do.asinh Logical. If \code{TRUE}, applies Arcsinh transformation using \code{do_asinh_local}.
#' @param do.asym Logical. If \code{TRUE}, applies an additional asymmetric root transformation to Arcsinh values.
#' @param asym_root Numeric. The root used for asymmetric transformation. Default is 2.
#' @param cofactor Numeric. Co-factor for Arcsinh transformation. Default is 500.
#' @param stim_to_exclude Character vector. List of stimulations to filter out from the final results.
#'
#' @return A \code{data.table} where each row is a single cell, containing metadata,
#'   gate positivity, and transformed intensity values.
#'
#' @import flowWorkspace flowCore tidyverse data.table
#' @export
extract_flow_exprs_data <- function(gs,
                                    output_nodes,
                                    parent_node,
                                    cytokine_nodes,
                                    pData_cols = c("BATCH", "PTID", "STIM", "VISITNO", "Run Num", "Replicate", "SAMP_ORD"),
                                    do.comp = FALSE,
                                    do.biexp = FALSE,
                                    do.asinh = TRUE,
                                    do.asym = TRUE,
                                    asym_root = 2,
                                    cofactor = 500,
                                    stim_to_exclude = NULL) {

  # 1. Pre-computation Checks
  # ---------------------------------------------------------------------------
  if (!do.asinh && do.asym) {
    stop("Arcsinh transformation must be TRUE to use Asymmetric transformation.")
  }

  # 2. Main Extraction Loop (Per Sample)
  # ---------------------------------------------------------------------------
  exprs.tmp <- flowWorkspace::lapply(gs, function(x) {

    # --- Metadata & Marker Mapping ---
    # Standardizes marker names and handles common naming discrepancies
    mark_names <- flowWorkspace::markernames(x)
    annotation <- data.frame(
      markername = names(mark_names),
      colname = sapply(mark_names, function(m) str_split(m, " ")[[1]][1]),
      row.names = names(mark_names)
    ) %>%
      dplyr::mutate(colname = str_replace_all(colname, "/", "_")) %>%
      dplyr::mutate(colname = case_when(
        colname == "Integrin" ~ "Integrin-B7",
        colname %in% c("Granzyme", "GzB") ~ "Granzyme-B",
        colname == "PD1" ~ "PD-1",
        colname == "TNFa" ~ "TNF",
        .default = colname
      ))

    # --- Data Extraction: Compensated & Biexponential ---
    # Get compensated data
    if(do.comp == TRUE){
      comp.FI <- flowCore::exprs(flowWorkspace::gh_pop_get_data(x, inverse.transform = TRUE))
      comp.FI <- comp.FI[, intersect(colnames(comp.FI), annotation$markername)]
      colnames(comp.FI) <- paste("comp", annotation[colnames(comp.FI), "colname"], sep = "_")
    }else{
      comp.FI <- NULL
    }
    # Get FlowJo biexp
    if(do.biexp == TRUE){
      biexp.FI <- flowCore::exprs(flowWorkspace::gh_pop_get_data(x, inverse.transform = FALSE))
      biexp.FI <- biexp.FI[, intersect(colnames(biexp.FI), annotation$markername)]
      colnames(biexp.FI) <- paste("biexp", annotation[colnames(biexp.FI), "colname"], sep = "_")
    }else{
      biexp.FI <- NULL
    }
    # Get arcsinh trans. data
    comp.FI.tmp <- flowCore::exprs(flowWorkspace::gh_pop_get_data(x, inverse.transform = TRUE))
    comp.FI.tmp <- comp.FI.tmp[, intersect(colnames(comp.FI.tmp), annotation$markername)]
    colnames(comp.FI.tmp) <- paste("comp", annotation[colnames(comp.FI.tmp), "colname"], sep = "_")
    markers <- colnames(comp.FI.tmp)
    asinh.FI <- ICSR::do_asinh_local(dat = comp.FI.tmp %>% as.data.table(), use.cols = markers, cofactor = cofactor)
    markers <- colnames(asinh.FI)[str_detect(string = colnames(asinh.FI), pattern = "asinh")]
    asinh.FI <- asinh.FI %>%
      dplyr::select(dplyr::all_of(markers))
    colnames(asinh.FI) <- str_replace_all(string = colnames(asinh.FI), pattern = "comp_", replacement = "asinh_")
    colnames(asinh.FI) <- str_remove(string = colnames(asinh.FI), pattern = "_asinh")
    if(do.asinh == FALSE){
      asinh.FI <- NULL
    }
    # Get arcsinh+asym trans. data
    asym_root_2_lo <- function(x, a = 2) {
      x <- case_when(x < (a - 1) ~ (a - abs(x-a)^(1/2)),
                     TRUE ~ x)
    }
    if(do.asym == TRUE){
      asinh.asym.FI <- apply(X = asinh.FI, MARGIN = 2, FUN = function(x) asym_root_2_lo(x = x, a = asym_root))
      colnames(asinh.asym.FI) <- str_replace(string = colnames(asinh.asym.FI), pattern = "asinh", replacement = "asinh_asym")
    }else{
      asinh.asym.FI <- NULL
    }

    # --- Population Positivity (Boolean Indices) ---
    # Retrieve true/false indices for each specified gate
    marker_response <- NULL
    # Using a retry loop for robust extraction from GatingSet
    while(is.null(marker_response)) {
      marker_response <- tryCatch({
        res <- lapply(output_nodes, function(mrkr) flowWorkspace::gh_pop_get_indices(x, mrkr))
        names(res) <- output_nodes
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
    dt.res <- dplyr::bind_cols(comp.FI, biexp.FI, asinh.FI, asinh.asym.FI, marker_response) %>%
      dplyr::bind_cols(pd_subset[rep(1, flowWorkspace::gh_pop_get_stats(x, "root")$count), ])

    # --- Final Filtering ---
    # Only keep cells belonging to the parent gate and expressing at least 1 cytokine
    dt.output <- dt.res %>%
      dplyr::filter(get(parent_node) == TRUE) %>%
      dplyr::mutate(NSUB = n())
    # Count how many cytokines are positive per cell
    dt.output$CYTNUM <- rowSums(dt.output[, cytokine_nodes, with = FALSE] == TRUE)
    # Filter for cytokine positive cells and select final columns
    final_cols <- c("FCS", standardized_pd_names, "NSUB", "CYTNUM", output_nodes,
                    grep("comp|biexp|asinh", colnames(dt.output), value = TRUE))
    dt.output <- dt.output %>%
      dplyr::filter(CYTNUM >= 1) %>%
      dplyr::select(dplyr::all_of(intersect(final_cols, colnames(dt.output))))
  })

  # 3. Final Aggregation and Filtering
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

