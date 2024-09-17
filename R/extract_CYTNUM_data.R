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


  ##-- Extraction
  exprs.tmp <- lapply(gs, function(x)
  {
    #- Require
    require(flowWorkspace)
    require(tidyverse)
    require(data.table)

    cat("\t", pData(x) %>% rownames(), "\n")

    #- Get boolean positivity call for markers
    options(warn = 0)
    marker_response <- try(lapply(c(parent_node, cytokine_nodes), function(mrkr){gh_pop_get_indices(x, mrkr)}))
    while(class(marker_response) == "try-error"){
      marker_response <- try(lapply(c(parent_node, cytokine_nodes), function(mrkr){gh_pop_get_indices(x, mrkr)}))
    }
    names(marker_response) <- c(parent_node, cytokine_nodes)
    marker_response <- bind_rows(marker_response)

    #- data.table()
    dt.res <- marker_response %>%
      data.table() %>%
      mutate(FCS = rownames(pData(x))) %>%
      mutate(BATCH = pData(x)$BATCH) %>%
      mutate(PTID = pData(x)$PTID) %>%
      mutate(STIM = pData(x)$STIM) %>%
      mutate(VISITNO = pData(x)$VISITNO) %>%
      mutate(RUNNUM = pData(x)$`Run Num`) %>%
      mutate(SAMP_ORD = pData(x)$SAMP_ORD) %>%
      mutate(REPLICATE = pData(x)$Replicate)
    dt.res$CYTNUM <- apply(dt.res[, ..cytokine_nodes], 1, function(x) sum(x == TRUE))

    #- Output
    dt.output <- dt.res %>%
      filter(get({{parent_node}}) == TRUE) %>%
      group_by(FCS, BATCH, PTID, STIM, VISITNO, RUNNUM, REPLICATE, SAMP_ORD) %>%
      mutate(NSUB = n()) %>%
      mutate(boolean_CYTNUM = factor(x = ifelse(CYTNUM >= 1, TRUE, FALSE), levels = c(TRUE, FALSE))) %>%
      group_by(BATCH, PTID, STIM, VISITNO, RUNNUM, REPLICATE, SAMP_ORD, NSUB, boolean_CYTNUM, .drop = FALSE) %>%
      summarize(CYTNUM = n())
    dt.output %>% return()
  })


  ##-- Output
  bind_rows(exprs.tmp) %>% return()
}

