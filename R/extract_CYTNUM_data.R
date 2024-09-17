#'
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Function
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' extract_CYTNUM_data
#'
#' @param gs GatingSet Object
#' @param parent_node Parent Flowjo gate
#' @param cytokine_nodes Cytokine Flowjo gates
#'
#' @return data.table
#' @export
#'
#' @examples
extract_CYTNUM_data <- function(gs,
                                parent_node,
                                cytokine_nodes)
{
  ##-- Require
  require(flowWorkspace)
  require(tidyverse)
  require(data.table)


  ##-- Extraction
  exprs.tmp <- lapply(gs, function(x)
  {
    cat("\t", flowWorkspace::pData(x) %>% rownames(), "\n")

    #- Get boolean positivity call for markers
    options(warn = 0)
    marker_response <- try(lapply(c(parent_node, cytokine_nodes), function(mrkr){flowWorkspace::gh_pop_get_indices(x, mrkr)}))
    while(class(marker_response) == "try-error"){
      marker_response <- try(lapply(c(parent_node, cytokine_nodes), function(mrkr){flowWorkspace::gh_pop_get_indices(x, mrkr)}))
    }
    names(marker_response) <- c(parent_node, cytokine_nodes)
    marker_response <- dplyr::bind_rows(marker_response)

    #- data.table()
    dt.res <- marker_response %>%
      data.table() %>%
      dplyr::mutate(FCS = rownames(flowWorkspace::pData(x))) %>%
      dplyr::mutate(BATCH = flowWorkspace::pData(x)$BATCH) %>%
      dplyr::mutate(PTID = flowWorkspace::pData(x)$PTID) %>%
      dplyr::mutate(STIM = flowWorkspace::pData(x)$STIM) %>%
      dplyr::mutate(VISITNO = flowWorkspace::pData(x)$VISITNO) %>%
      dplyr::mutate(RUNNUM = flowWorkspace::pData(x)$`Run Num`) %>%
      dplyr::mutate(SAMP_ORD = flowWorkspace::pData(x)$SAMP_ORD) %>%
      dplyr::mutate(REPLICATE = flowWorkspace::pData(x)$Replicate)
    dt.res$CYTNUM <- apply(dt.res[, ..cytokine_nodes], 1, function(x) sum(x == TRUE))

    #- Output
    dt.output <- dt.res %>%
      dplyr::filter(get({{parent_node}}) == TRUE) %>%
      dplyr::group_by(FCS, BATCH, PTID, STIM, VISITNO, RUNNUM, REPLICATE, SAMP_ORD) %>%
      dplyr::mutate(NSUB = n()) %>%
      dplyr::mutate(boolean_CYTNUM = factor(x = ifelse(CYTNUM >= 1, TRUE, FALSE), levels = c(TRUE, FALSE))) %>%
      dplyr::group_by(BATCH, PTID, STIM, VISITNO, RUNNUM, REPLICATE, SAMP_ORD, NSUB, boolean_CYTNUM, .drop = FALSE) %>%
      dplyr::summarize(CYTNUM = n())
    dt.output %>% return()
  })


  ##-- Output
  dplyr::bind_rows(exprs.tmp) %>% return()
}

