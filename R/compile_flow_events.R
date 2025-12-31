#' Compile Flow Cytometry Events and Polyfunctionality Statistics
#'
#' This wrapper function serves as the primary entry point for data extraction in the \code{ICSR}
#' pipeline. It simultaneously extracts detailed single-cell expression data and
#' summarized cytokine polyfunctionality counts by calling \code{extract_flow_exprs_data}
#' and \code{extract_CYTNUM_data}.
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
#' @return A \code{list} containing two \code{data.table} elements:
#' \itemize{
#'   \item \code{exprs}: Single-cell level intensities and gate positivity.
#'   \item \code{cytnum}: Summarized counts of cytokine-positive vs. negative cells per sample.
#' }
#'
#' @export
compile_flow_events <- function(gs,
                                output_nodes,
                                parent_node,
                                cytokine_nodes,
                                pData_cols = c("BATCH", "PTID", "STIM", "VISITNO", "Run Num", "Collection Num", "Replicate", "SAMP_ORD"),
                                do.comp = FALSE,
                                do.biexp = FALSE,
                                do.asinh = TRUE,
                                do.asym = TRUE,
                                asym_root = 2,
                                cofactor = 500,
                                stim_to_exclude = NULL)
{
  # 1. User Notifications
  # ---------------------------------------------------------------------------
  if (!is.null(stim_to_exclude)) {
    message("Filtering data: The following stimulations will be excluded: ",
            paste(stim_to_exclude, collapse = ", "))
  }

  # 2. Extract Single-Cell Expression Data
  # ---------------------------------------------------------------------------
  # This returns intensities (Comp, Biexp, Asinh) and Boolean calls
  message("Step 1/2: Extracting single-cell expression data...")
  dt.exprs <- ICSR::extract_flow_exprs_data(gs = gs,
                                            output_nodes = output_nodes,
                                            parent_node = parent_node,
                                            cytokine_nodes = cytokine_nodes,
                                            pData_cols = c("BATCH", "PTID", "STIM", "VISITNO", "Run Num", "Collection Num", "Replicate", "SAMP_ORD"),
                                            do.comp = do.comp,
                                            do.biexp = do.biexp,
                                            do.asinh = do.asinh,
                                            do.asym = do.asym,
                                            asym_root = asym_root,
                                            cofactor = cofactor,
                                            stim_to_exclude = stim_to_exclude)

  # 3. Extract Polyfunctionality (CYTNUM) Data
  # ---------------------------------------------------------------------------
  # This returns summarized counts of cytokine-expressing cells
  message("Step 2/2: Extracting polyfunctionality summary statistics...")
  dt.cytnum <- ICSR::extract_CYTNUM_data(gs = gs,
                                         parent_node = parent_node,
                                         cytokine_nodes = cytokine_nodes,
                                         pData_cols = c("BATCH", "PTID", "STIM", "VISITNO", "Run Num", "Collection Num", "Replicate", "SAMP_ORD"),
                                         stim_to_exclude = stim_to_exclude)

  # 4. Final Output Construction
  # ---------------------------------------------------------------------------
  #
  return(list(
    "exprs" = dt.exprs,
    "cytnum" = dt.cytnum
  ))
}

