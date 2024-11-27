#'
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Function
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Filtering
#'
#' @param dt input data.table. Including pData and fluorescence intensities (FI).
#' @param remove_dup TRUE or FALSE. If TRUE, remove duplicates (based on RUNNUM).
#' @param dt.MTL_dup Master thaw list to remove duplicates (based on RUNNUM). Need BATCH, SAMP_ORD, PTID, VISITNO and RUNNUM.
#' @param remove_participants TRUE or FALSE. If TRUE, remove participants, e.g. control participants.
#' @param participants_list Vector of removed participants.
#' @param remove_nsub TRUE or FALSE. If TRUE, remove samples based on total number of cells. Ex: 10,000 CD4+ T.
#' @param cutoff_nsub Cut-off for `remove_nsub`. Ex: 10,000 CD4+ T.
#' @param flag_unreliable TRUE or FALSE. Flag unreliable samples. If TRUE, need unreliable lists.
#' @param unreliable_list_all Vector of unreliable samples. Unreliable sample level must be BATCH_PTID_VISITNO_RUNNUM. Use with caution.
#' @param unreliable_list_by_stim Vector of unreliable samples. Unreliable sample level must be BATCH_PTID_VISITNO_RUNNUM_REPLICATE_STIM. Use with caution.
#'
#' @return data.table with pData, FlowJo positivity (positivity for each marker using the marker gates determined by the lab (FlowJo gating)), and FI after filtering.
#'
#' @export
#'
#' @examples
#'
do_filtering <- function(dt,
                         remove_dup = FALSE,
                         dt.MTL_dup = NULL,
                         remove_participants = FALSE,
                         participants_list = NULL,
                         remove_nsub = FALSE,
                         cutoff_nsub = NULL,
                         flag_unreliable = FALSE,
                         unreliable_list_all = NULL,
                         unreliable_list_by_stim = NULL) {
  #- Require
  require(tidyverse)

  #- Checks
  if(remove_dup == TRUE & is.null(dt.MTL_dup)){
    stop("do_filtering needs Master Thaw List to remove duplicates.")
  }

  #- Summary - duplicates
  # Needs MTL file
  if(!is.null(dt.MTL_dup)){
    batch.df <- dt.MTL_dup %>%
      dplyr::select(BATCH, SAMP_ORD, PTID, VISITNO, RUNNUM)
    batch.df <- batch.df %>%
      dplyr::group_by(PTID, VISITNO) %>%
      dplyr::mutate(dup = ifelse(RUNNUM == max(RUNNUM), FALSE, TRUE))
  }

  #- Remove duplicates
  if(remove_dup == TRUE){
    dups <- batch.df %>%
      dplyr::filter(dup == TRUE) %>%
      dplyr::mutate(SAMPLE = paste(BATCH, SAMP_ORD, PTID, VISITNO, RUNNUM)) %>%
      dplyr::pull(SAMPLE)
    dt <- dt %>%
      dplyr::filter(!(paste(BATCH, SAMP_ORD, PTID, VISITNO, RUNNUM) %in% dups))
  }

  #- Remove participants
  if(remove_participants == TRUE){
    dt <- dt %>%
      dplyr::filter(!(PTID %in% participants_list))
  }

  #- Remove based on nsub
  if(remove_nsub == TRUE){
    dt <- dt %>%
      dplyr::filter(NSUB > cutoff_nsub)
  }

  #- Flag unreliable
  if(flag_unreliable == TRUE){
    dt <- dt %>%
      dplyr::mutate(reliable_flag = dplyr::case_when(paste(BATCH, PTID, VISITNO, RUNNUM, sep = "_") %in% unreliable.s_1 ~ "0",
                                                     paste(BATCH, PTID, VISITNO, RUNNUM, REPLICATE, STIM, sep = "_") %in% unreliable.s_2 ~ "0",
                                                     .default = "1"))
  }

  #- Output
  dt %>%
    return()
}

