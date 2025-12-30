#' Create and Annotate GatingSet from Raw Data
#'
#' This function parses a FlowJo XML workspace, associates it with raw FCS files,
#' and merges experimental metadata (thaw list/batch data) into the GatingSet's
#' phenotypic data (pData).
#'
#' @param assayid Character. The unique ID for the assay (must match a value in `batchData$Batch`).
#' @param xml_path Character. Path to the FlowJo XML workspace file.
#' @param fcs_path Character. Path to the directory containing the FCS files.
#' @param sample_name Character. The name of the sample group/node in the XML to load. Default is "Samples".
#' @param xml_keywords Character vector. Keywords to extract from the XML workspace.
#'
#' @return Saves a GatingSet object to `data-raw/tmpdata/[assayid]` and returns the path invisibly.
#'
#' @details
#' The function relies on a global or parent-environment object `batchData` which
#' must contain `Batch` and `SAMP_ORD` columns. The resulting GatingSet is saved
#' using `flowWorkspace::save_gs`.
#'
#' @examples
raw_to_GatingSet <- function(assayid,
                             xml_path,
                             fcs_path,
                             sample_name = "Samples",
                             xml_keywords = c("$FIL", "$TOT", "EXPERIMENT NAME", "Sample Order", "Stim", "Replicate"))
{
  message("Processing BATCH: ", assayid)

  # 1. Parse FlowJo XML and convert to GatingSet
  # ---------------------------------------------------------------------------
  workspace <- CytoML::open_flowjo_xml(file = xml_path)
  G <- CytoML::flowjo_to_gatingset(workspace,
                                   name = sample_name,
                                   keywords = xml_keywords,
                                   path = fcs_path,
                                   additional.sampleID = TRUE)

  # 2. Extract and format Phenotypic Data (pData)
  # ---------------------------------------------------------------------------
  # Filter the master thaw list for the current batch
  batchdf <- batchData[which(batchData$Batch == assayid), ]
  pd <- flowWorkspace::pData(G)
  pd$BATCH <- unique(batchdf$Batch)
  # Ensure SAMP_ORD is numeric to allow for successful merging with batchData
  # Note: Case-insensitive search for the 'Sample Order' keyword extracted from XML
  pd$SAMP_ORD <- as.numeric(as.character(pd[, which(tolower(names(pd)) %in% "sample order")]))
  pd$STIM <- pd[, which(tolower(names(pd)) %in% "stim")]
  pd$roworder <- 1:nrow(pd)
  pd$sample_name <- rownames(pd)

  # 3. Merge Metadata
  # ---------------------------------------------------------------------------
  # Join the GatingSet pData with the Batch/Thaw information
  pd <- merge(x = pd[, c("name", "SAMP_ORD", "BATCH", "STIM", "Replicate", "EXPERIMENT NAME", "roworder", "sample_name")],
              y = batchdf,
              by.x = c("SAMP_ORD", "BATCH"),
              by.y = c("SAMP_ORD", "Batch"),
              all.x = TRUE)
  # Restore original row names and sort order after merge
  rownames(pd) <- pd[["sample_name"]]
  pd <- pd[order(pd$roworder), ]
  flowWorkspace::pData(G) <- pd

  # 4. Save and Export
  # ---------------------------------------------------------------------------
  # Ensure the directory /data-raw/tmpdata/ exists before calling this function
  save_path <- here::here("data-raw", "tmpdata", assayid)
  flowWorkspace::save_gs(G, path = save_path, overwrite = TRUE)
  message("GatingSet successfully saved to: ", save_path)
}

