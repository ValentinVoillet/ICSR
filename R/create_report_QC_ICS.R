#'
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Function
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' ICS QC Report
#'
#' This function creates a data QC report.
#' @param dt.exprs input data.table. Including FI.
#' @param dt.cytnum input data.table. Including # of cytokines per sample.
#' @param output_format output format in \link{render}. Default is \code{html_document(toc = TRUE, toc_depth = 6, theme = "yeti")}.
#' @param output_file output file name in \link{render}. Default is "report.html".
#' @param output_dir output directory for report in \link{render}. Default is user's current directory.
#' @param report_title report title. Default is "ICS QC Report".
#' @param report_author report authors. Default is NULL.
#'
#' @return
#'
#' @export
#'
#' @examples
#'
create_report_QC_ICS <- function(dt.exprs,
                                 dt.cytnum,
                                 output_format = html_document(toc = TRUE, toc_depth = 6, theme = "yeti"),
                                 output_file = "ICS_QC_report.html",
                                 output_dir = getwd(),
                                 report_title = "ICS QC Report",
                                 report_author = NULL) {
  require("rmarkdown")
  require("data.table")

  #- Check if input is data.table
  if (!data.table::is.data.table(dt.exprs)) dt.exprs <- data.table::data.table(dt.exprs)
  if (!data.table::is.data.table(dt.cytnum)) dt.cytnum <- data.table::data.table(dt.cytnum)

  #- Get directory of report markdown template
  # report_dir <- ("rmd_template/QC_report.rmd")
  report_dir <- system.file("rmd_template/QC_report.rmd", package = "ICSR")

  #- Render report into html
  suppressWarnings(rmarkdown::render(
    input = report_dir,
    output_format = output_format,
    output_file = output_file,
    output_dir = output_dir,
    intermediates_dir = output_dir,
    params = list(dt.exprs = dt.exprs, dt.cytnum = dt.cytnum, set_title = report_title, set_author = report_author)))

  #- Open report
  report_path <- path.expand(file.path(output_dir, output_file))
  browseURL(report_path)
}

