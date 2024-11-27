#'
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Function
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' create_gating_set
#'
#' @param assayid `numeric` or `character`. batch ID: used for merging Master thaw list and `pData()`.
#' @param xml_path xml path (FlowJo .xml).
#' @param fcs_path fcs path.
#' @param batchData Master thaw list. Need BATCH, SAMP_ORD, PTID, VISITNO and RUNNUM.
#' @param sample_name `numeric` or `character`. The name or index of the group of samples to be imported.
#' @param xml_keywords `character` vector specifying the keywords to be extracted as pData of GatingSet.
#' @param output_path Output path for the new GatingSet object.
#'
#' @return GatingSet R obj.
#'
#' @export
#'
#' @examples
create_gating_set <- function(assayid,
                              xml_path,
                              fcs_path,
                              batchData,
                              sample_name = "Samples",
                              xml_keywords = c("$FIL", "$TOT", "EXPERIMENT NAME", "Sample Order", "Stim", "Replicate"),
                              output_path = NULL) {

  #- Require
  require(CytoML)
  require(flowWorkspace)

  #- GatingSet
  workspace <- CytoML::open_flowjo_xml(file = xml_path)
  G <- CytoML::flowjo_to_gatingset(workspace,
                                   name = sample_name,
                                   keywords = xml_keywords,
                                   path = fcs_path,
                                   additional.sampleID = TRUE)

  #- Add annotation data to gating set
  batchdf <- batchData[which(batchData$BATCH == assayid), ]
  pd <- flowWorkspace::pData(G)
  pd$BATCH <- unique(batchdf$BATCH)
  # Find Sample Order column and convert to numeric. SAMP_ORD is key for merge with batch data
  pd$SAMP_ORD <- as.numeric(as.character(pd[, which(tolower(names(pd)) %in% "sample order")]))
  pd$STIM <- pd[, which(tolower(names(pd)) %in% "stim")]
  pd$roworder <- 1:nrow(pd)
  pd$sample_name <- rownames(pd)

  #- Add batchData information to pData, then apply to GatingSet.
  pd <- merge(x = pd[,c("name", "SAMP_ORD", "BATCH", "STIM", "Replicate", "EXPERIMENT NAME", "roworder", "sample_name")],
              y = batchdf,
              by.x = c("SAMP_ORD", "BATCH"),
              by.y = c("SAMP_ORD", "BATCH"),
              all.x = TRUE)
  rownames(pd) <- pd[["sample_name"]]
  pd <- pd[order(pd$roworder), ]
  flowWorkspace::pData(G) <- pd

  #- Output
  path <- paste0(output_path, "/", assayid)
  flowWorkspace::save_gs(G, path = path, overwrite = TRUE)
}

