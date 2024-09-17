#' @include do_asinh_local.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Function
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' extract_flow_exprs_data
#'
#' @param gs GatingSet Object
#' @param output_nodes Flowjo gates
#' @param parent_node Parent Flowjo gate
#' @param cytokine_nodes Cytokine Flowjo gates
#' @param do.comp Get compensated values
#' @param do.biexp Get biexp transformed values using FlowJo transformation parameters
#' @param do.asinh Get arcsinh transformed values
#' @param cofactor Cofactor used for arcsinh transformation
#'
#' @return data.table
#' @export
#'
#' @examples
extract_flow_exprs_data <- function(gs,
                               output_nodes,
                               parent_node,
                               cytokine_nodes,
                               do.comp = TRUE,
                               do.biexp = TRUE,
                               do.asinh = TRUE,
                               cofactor = 500)
{
  ##-- Require
  require(flowWorkspace)
  require(tidyverse)


  ##-- Extraction
  exprs.tmp <- flowWorkspace::lapply(gs, function(x)
  {
    cat("\t", pData(x) %>% rownames(), "\n")

    #- Annotation
    annotation <- data.frame(markername = flowWorkspace::markernames(x) %>% names(),
                             colname = sapply(X = flowWorkspace::markernames(x), FUN = function(x) str_split(string = x, pattern = " ")[[1]][1]),
                             row.names = flowWorkspace::markernames(x) %>% names()) %>%
      dplyr::mutate(colname = str_replace_all(string = colname, pattern = "/", replacement = "_")) %>%
      dplyr::mutate(colname = case_when(colname == "Integrin" ~ "Integrin-B7",
                                        colname == "Granzyme" ~ "Granzyme-B",
                                        colname == "GzB" ~ "Granzyme-B",
                                        colname == "PD1" ~ "PD-1",
                                        colname == "TNFa" ~ "TNF",
                                        .default = colname))

    #- Get compensated data
    if(do.comp == TRUE){
      comp.FI <- flowCore::exprs(flowWorkspace::gh_pop_get_data(x, inverse.transform = TRUE))
      comp.FI <- comp.FI[, intersect(colnames(comp.FI), annotation$markername)]
      colnames(comp.FI) <- paste("comp", annotation[colnames(comp.FI), "colname"], sep = "_")
    }else{
      comp.FI <- NULL
    }

    #- Get FlowJo biexp
    if(do.biexp == TRUE){
      biexp.FI <- flowCore::exprs(flowWorkspace::gh_pop_get_data(x, inverse.transform = FALSE))
      biexp.FI <- biexp.FI[, intersect(colnames(biexp.FI), annotation$markername)]
      colnames(biexp.FI) <- paste("biexp", annotation[colnames(biexp.FI), "colname"], sep = "_")
    }else{
      biexp.FI <- NULL
    }

    #- Get arcsinh trans. data
    if(do.asinh == TRUE){
      comp.FI <- flowCore::exprs(flowWorkspace::gh_pop_get_data(x, inverse.transform = TRUE))
      comp.FI <- comp.FI[, intersect(colnames(comp.FI), annotation$markername)]
      colnames(comp.FI) <- paste("comp", annotation[colnames(comp.FI), "colname"], sep = "_")
      markers <- colnames(comp.FI)
      asinh.FI <- ICSR::do_asinh_local(dat = comp.FI %>% as.data.table(),
                                       use.cols = markers,
                                       cofactor = cofactor)
      markers <- colnames(asinh.FI)[str_detect(string = colnames(asinh.FI), pattern = "asinh")]
      asinh.FI <- asinh.FI[, ..markers]
      colnames(asinh.FI) <- str_replace_all(string = colnames(asinh.FI), pattern = "comp_", replacement = "asinh_")
      colnames(asinh.FI) <- str_remove(string = colnames(asinh.FI), pattern = "_asinh")
    }else{
      asinh.FI <- NULL
    }

    #- Get boolean positivity call for markers
    options(warn = 0)
    marker_response <- try(lapply(output_nodes, function(mrkr){flowWorkspace::gh_pop_get_indices(x, mrkr)}))
    while(class(marker_response) == "try-error"){
      marker_response <- try(lapply(output_nodes, function(mrkr){flowWorkspace::gh_pop_get_indices(x, mrkr)}))
    }
    names(marker_response) <- output_nodes
    marker_response <- dplyr::bind_rows(marker_response)

    #- data.table()
    dt.res <- dplyr::bind_cols(comp.FI, biexp.FI, asinh.FI, marker_response) %>%
      dplyr::mutate(FCS = rownames(pData(x))) %>%
      dplyr::mutate(BATCH = pData(x)$BATCH) %>%
      dplyr::mutate(PTID = pData(x)$PTID) %>%
      dplyr::mutate(STIM = pData(x)$STIM) %>%
      dplyr::mutate(VISITNO = pData(x)$VISITNO) %>%
      dplyr::mutate(RUNNUM = pData(x)$`Run Num`) %>%
      dplyr::mutate(REPLICATE = pData(x)$Replicate) %>%
      dplyr::mutate(SAMP_ORD = pData(x)$SAMP_ORD)

    #- Output
    dt.output <- dt.res %>%
      dplyr::filter(get({{parent_node}}) == TRUE) %>%
      dplyr::group_by(FCS, BATCH, PTID, STIM, VISITNO, RUNNUM, REPLICATE, SAMP_ORD) %>%
      dplyr::mutate(NSUB = n())
    dt.output$CYTNUM <- apply(dt.output[, cytokine_nodes], 1, function(x) sum(x == TRUE))
    cols <- c("FCS", "BATCH", "PTID", "STIM", "VISITNO", "RUNNUM", "REPLICATE", "SAMP_ORD", "NSUB", "CYTNUM",
              colnames(dt.output)[str_detect(string = colnames(dt.output), pattern = "Time")],
              colnames(dt.output)[str_detect(string = colnames(dt.output), pattern = "comp")],
              colnames(dt.output)[str_detect(string = colnames(dt.output), pattern = "biexp")],
              colnames(dt.output)[str_detect(string = colnames(dt.output), pattern = "asinh")])
    dt.output <- dt.output %>%
      dplyr::ungroup() %>%
      dplyr::filter(CYTNUM >= 1) %>%
      dplyr::select(all_of(cols))
    dt.output %>% return()
  })


  ##-- Output
  dplyr::bind_rows(exprs.tmp) %>% return()
}

