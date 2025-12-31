#' Generate an Automated ICS Quality Control Report
#'
#' This function takes the extracted expression and polyfunctionality data to
#' generate a comprehensive HTML QC report. The report includes visualizations for
#' marker distributions (ridge plots), cytokine positivity, and sample-level
#' statistics to guide data filtering.
#'
#' @param dt.exprs A \code{data.table} containing single-cell expression data
#'   (typically the \code{exprs} element from \code{compile_flow_events}).
#' @param dt.cytnum A \code{data.table} containing summarized population counts
#'   (typically the \code{cytnum} element from \code{compile_flow_events}).
#' @param cytokine_nodes Character vector. Names of cytokine gates to include
#'   in the polyfunctionality summary.
#' @param markers Character vector. List of markers to visualize using ridge plots.
#'   Should match column names in \code{dt.exprs}.
#' @param output_format The RMarkdown render format. Defaults to a
#'   \code{html_document} with a floating TOC and the 'yeti' theme.
#' @param output_file Character. The name of the resulting HTML file.
#'   Default is "ICS_QC_report.html".
#' @param output_dir Character. The directory where the report will be saved.
#'   Defaults to the current working directory.
#' @param report_title Character. The title displayed at the top of the report.
#' @param report_author Character. The name of the author/lab to display.
#'
#' @return Invisibly returns the path to the generated report and opens it in
#'   the default web browser.
#'
#' @import rmarkdown data.table
#' @export
#'
#' @examples
#' # create_report_QC_ICS(dt.exprs = output$exprs, dt.cytnum = output$cytnum)
create_report_QC_ICS <- function(dt.exprs,
                                 dt.cytnum,
                                 cytokine_nodes = NULL,
                                 markers = NULL,
                                 output_format = rmarkdown::html_document(toc = TRUE, toc_depth = 6, theme = "yeti"),
                                 output_file = "ICS_QC_report.html",
                                 output_dir = getwd(),
                                 report_title = "ICS QC Report",
                                 report_author = NULL) {
  # 1. Input Validation
  # ---------------------------------------------------------------------------
  # Ensure inputs are data.tables to support fast processing within the Rmd templateif(!data.table::is.data.table(dt.exprs)) dt.exprs <- data.table::data.table(dt.exprs)
  if(!data.table::is.data.table(dt.cytnum)) dt.cytnum <- data.table::data.table(dt.cytnum)

  # 2. Locate Template
  # ---------------------------------------------------------------------------
  # Finds the RMarkdown template bundled within the ICSR package
  report_dir <- system.file("rmd_template/QC_report.rmd", package = "ICSR")

  # 3. Render Report
  # ---------------------------------------------------------------------------
  # Pass data and parameters into the Rmd file as a list
  message("Generating HTML QC report...")
  suppressWarnings(rmarkdown::render(
    input = report_dir,
    output_format = output_format,
    output_file = output_file,
    output_dir = output_dir,
    intermediates_dir = output_dir,
    params = list(dt.exprs = dt.exprs,
                  dt.cytnum = dt.cytnum,
                  cytokine_nodes = cytokine_nodes,
                  markers = markers,
                  set_title = report_title,
                  set_author = report_author),
    quiet = TRUE # Set to FALSE if debugging template errors
  ))

  # 4. Finalization
  # ---------------------------------------------------------------------------
  # Build the absolute path and automatically open the report for review
  report_path <- path.expand(file.path(output_dir, output_file))
  if (file.exists(report_path)) {
    message("Report successfully generated: ", report_path)
    utils::browseURL(report_path)
  } else {
    warning("Report file was not found after rendering.")
  }
  browseURL(report_path)
}

