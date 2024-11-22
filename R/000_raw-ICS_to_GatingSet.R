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


library(here)
library(tidyverse)
library(data.table)
library(doMC)
library(XML)
library(flowWorkspace)
library(CytoML)
library(parallel)

#- GatingSet function
create_gating_set <- function(assayid,
                              xml_path,
                              fcs_path,
                              sample_name = "Samples",
                              xml_keywords = c("$FIL", "$TOT", "EXPERIMENT NAME", "Sample Order", "Stim", "Replicate")) {


  cat(assayid, "\n")

  # GatingSet
  workspace <- CytoML::open_flowjo_xml(file = xml_path)
  G <- CytoML::flowjo_to_gatingset(workspace,
                                   name = sample_name,
                                   keywords = xml_keywords,
                                   path = fcs_path,
                                   additional.sampleID = TRUE)

  # Add annotation data to gating set
  batchdf <- batchData[which(batchData$Batch == assayid), ]
  pd <- flowWorkspace::pData(G)
  pd$BATCH <- unique(batchdf$Batch)
  # Find Sample Order column and convert to numeric. SAMP_ORD is key for merge with batch data
  pd$SAMP_ORD <- as.numeric(as.character(pd[, which(tolower(names(pd)) %in% "sample order")]))
  pd$STIM <- pd[, which(tolower(names(pd)) %in% "stim")]
  pd$roworder <- 1:nrow(pd)
  pd$sample_name <- rownames(pd)

  # Add batchData information to pData, then apply to GatingSet.
  pd <- merge(x = pd[,c("name", "SAMP_ORD", "BATCH", "STIM", "Replicate", "EXPERIMENT NAME", "roworder", "sample_name")],
              y = batchdf,
              by.x = c("SAMP_ORD", "BATCH"),
              by.y = c("SAMP_ORD", "Batch"),
              all.x = TRUE)
  rownames(pd) <- pd[["sample_name"]]
  pd <- pd[order(pd$roworder), ]
  flowWorkspace::pData(G) <- pd

  # Output
  #flowWorkspace::save_gs(G, path = here::here("data-raw", "tmpdata", assayid), overwrite = TRUE)
  flowWorkspace::save_gs(G, path = here::here("data-raw", "version_2", "tmpdata", assayid), overwrite = TRUE)
}


