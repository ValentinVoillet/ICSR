#'
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Function
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' dotplot_annotation
#'
#' @param dt data.table Input sample.
#' @param flowjo_gates flowJo gates.
#' @param markers List of markers - MFI calculation.
#' @param cluster_col Cluster column.
#' @param level_marker marker order to follow.
#' @param order cluster order to follow.
#'
#' @return UMAP embeddings
#' @export
#'
#' @examples
dotplot_annotation <- function(dt, markers, flowjo_gates, cluster_col = "Leiden", level_marker = NULL, order = NULL){
  #- dt.ggplot_boolean
  dt.boolean <- dt.ggplot_boolean(dt, flowjo_gates, cluster_col = "Leiden")

  #- dt.ggplot_MFI
  dt.MFI <- dt.ggplot_MFI(dt, markers, cluster_col = "Leiden")

  #- Annotations
  x_axis_annotation <- x_axis(dt.boolean = dt.boolean, order = order)

  #- Processing - dot-plot
  dt.boolean$variable <- str_remove_all(string = dt.boolean$variable, pattern = "[+]")
  dt.MFI$variable <- str_remove_all(string = dt.MFI$variable, pattern = "asinh_asym_|asinh|biexp_|comp_")
  dt.ggplot_merged <- merge(x = dt.boolean %>%
                              ungroup() %>%
                              select(cluster, variable, marker_positivity),
                            y = dt.MFI %>%
                              select(cluster, variable, MFI),
                            all.x = TRUE,
                            all.y = TRUE)
  dt.ggplot_merged <- dt.ggplot_merged %>%
    dplyr::mutate(facet = case_when(is.na(MFI) ~ "A",
                                    is.na(marker_positivity) ~ "C",
                                    .default = "B"))
  dt.ggplot_merged <- dt.ggplot_merged %>%
    mutate(marker_positivity = case_when(is.na(marker_positivity) ~ 50, .default = marker_positivity))
  if(is.null(order) == FALSE){
    dt.ggplot_merged <- dt.ggplot_merged %>%
      dplyr::mutate(cluster = factor(x = cluster, levels = order))
  }
  if(is.null(level_marker) == FALSE){
    dt.ggplot_merged <- dt.ggplot_merged %>%
      dplyr::mutate(variable = factor(x = variable, levels = level_marker))
  }

  #- Dot-plot
  dt.ggplot_merged %>%
    dplyr::mutate(variable = fct_rev(f = variable)) %>%
    ggplot(aes(x = cluster, y = variable)) +
      geom_point(aes(size = marker_positivity, color = MFI)) +
      geom_text(aes(label = round(marker_positivity, 0)),
                angle = 0, size = 5, color = "white",
                data = dt.ggplot_merged %>% dplyr::filter(marker_positivity > 20 & facet != "C")) +
      facet_grid(rows = vars(facet), scales = "free_y", space = "free_y") +
      scale_color_viridis_c(option = "D", limits = c(-1, 8)) +
      scale_x_discrete(labels = x_axis_annotation) +
      scale_size(range = c(1, 20), limits = c(0, 100)) +
      labs(color = "Median arcsinh(x/500) value", size = "Percent positive, gated") +
      theme_bw() +
      theme(axis.text.x = element_text(size = 14, angle = 60, hjust = 1, vjust = 1),
            axis.text.y = element_text(size = 20),
            axis.title = element_blank(),
            panel.grid = element_blank(),
            strip.background = element_blank(),
            strip.text = element_blank(),
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 18),
            legend.position = "bottom",
            legend.box = "vertical") -> plot
  plot %>%
    return()
}

