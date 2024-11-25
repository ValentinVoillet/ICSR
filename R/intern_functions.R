#'
#'
NULL

#' dt.ggplot_boolean
#'
#' @param dt data.table Input sample.
#' @param flowjo_gates flowJo gates.
#' @param cluster_col Cluster column.
#'
#' @return data.table with marker positivity, gated
#' @export
#'
#' @examples
dt.ggplot_boolean <- function(dt, flowjo_gates, cluster_col = "Leiden"){
  #- n's
  dt.summary.tmp <- dt %>%
    dplyr::group_by(get({{cluster_col}})) %>%
    dplyr::summarise(n = n())
  n <- dt.summary.tmp$n
  names(n) <- dt.summary.tmp[, 1] %>% unlist()

  #- % positive, gated & output
  dt.tmp <- dt %>%
    dplyr::select(dplyr::all_of(c(cluster_col, flowjo_gates))) %>%
    reshape2::melt(id.vars = cluster_col) %>%
    dplyr::mutate(variable = factor(x = variable, levels = flowjo_gates)) %>%
    dplyr::group_by(get({{cluster_col}}), value, variable, .drop = FALSE) %>%
    dplyr::summarize(n = n()) %>%
    dplyr::ungroup()
  colnames(dt.tmp)[1] <- "cluster"
  dt.tmp$n.total <- plyr::mapvalues(x = dt.tmp$cluster, from = names(n), to = n, warn_missing = FALSE)
  dt.tmp %>%
    dplyr::mutate(n.total = as.numeric(as.character(n.total))) %>%
    dplyr::filter(value == TRUE) %>%
    dplyr::mutate(marker_positivity = 100 * n / n.total) %>%
    data.table::data.table() %>%
    return()
}

#' dt.ggplot_MFI
#'
#' @param dt data.table Input sample.
#' @param markers List of markers - MFI calculation.
#' @param cluster_col Cluster column.
#'
#' @return data.table with MFI
#' @export
#'
#' @examples
dt.ggplot_MFI <- function(dt, markers, cluster_col = "Leiden"){
  #- MFI
  dt.tmp <- dt %>%
    dplyr::ungroup() %>%
    dplyr::select(cluster_col, markers) %>%
    reshape2::melt(id.vars = cluster_col) %>%
    dplyr::group_by(Leiden, variable) %>%
    dplyr::summarize(MFI = median(value))
  colnames(dt.tmp)[1] <- "cluster"

  #- Outpu
  dt.tmp %>%
    data.table::data.table() %>%
    return()
}

#' x_axis
#'
#' @param dt.boolean data.table boolean input.
#' @param order cluster order to follow.
#'
#' @return Vector for x axis.
#' @export
#'
#' @examples
x_axis <- function(dt.boolean, order = NULL) {
  #- 80%, 60% and 40%
  pos80lbl <- dt.boolean %>%
    dplyr::filter(marker_positivity >= 80) %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(annot = paste0("80% (", paste0(variable, collapse = ""), ")"))
  pos60lbl <- dt.boolean %>%
    dplyr::filter(marker_positivity < 80) %>%
    dplyr::filter(marker_positivity >= 60) %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(annot = paste0("60% (", paste0(variable, collapse = ""), ")"))
  pos40lbl <- dt.boolean %>%
    dplyr::filter(marker_positivity < 60) %>%
    dplyr::filter(marker_positivity >= 40) %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(annot = paste0("40% (", paste0(variable, collapse = ""), ")"))

  #- Output
  dt.annotation <- merge(x = pos80lbl, y = pos60lbl, by = "cluster", all = TRUE)
  dt.annotation <- merge(x = dt.annotation, y = pos40lbl, by = "cluster", all = TRUE)
  annotation <- dt.annotation %>%
    mutate(axis.text = paste(cluster, annot.x, annot.y, annot, sep = "\n")) %>%
    mutate(axis.text = str_remove(string = axis.text, pattern = "80% [(][])]\n")) %>%
    mutate(axis.text = str_remove(string = axis.text, pattern = "60% [(][])]\n")) %>%
    mutate(axis.text = str_remove(string = axis.text, pattern = "40% [(][])]")) %>%
    mutate(axis.text = str_remove(string = axis.text, pattern = "\nNA")) %>%
    mutate(cluster = factor(x = cluster, levels = order)) %>%
    arrange(cluster) %>%
    pull(axis.text)
  annotation %>% return()
}

