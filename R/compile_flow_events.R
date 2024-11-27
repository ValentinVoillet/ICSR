#' @include do_asinh_local.R
#' @include extract_flow_exprs_data.R
#' @include extract_CYTNUM_data.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Function
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' compile_flow_events
#'
#' @param gs GatingSet Object.
#' @param output_nodes Flowjo gates.
#' @param parent_node Parent Flowjo gate.
#' @param cytokine_nodes Cytokine Flowjo gates.
#' @param do.comp TRUE or FALSE. If TRUE, get compensated values.
#' @param do.biexp TRUE or FALSE. If TRUE, get biexp transformed values using FlowJo transformation parameters.
#' @param do.asinh TRUE or FALSE. If TRUE, get arcsinh transformed values.
#' @param do.asinh TRUE or FALSE. If TRUE, get arcsinh + asym transformed values.
#' @param asym_root Root used for asym transformation. By default 2.
#' @param cofactor Co-factor used for arcsinh transformation.
#' @param stim_to_exclude Stimulation(s) to be removed. By default NULL.
#'
#' @return list with two data.table elements 'exprs' and 'cytnum'
#'
#' @export
#'
#' @examples
compile_flow_events <- function(gs,
                                output_nodes,
                                parent_node,
                                cytokine_nodes,
                                do.comp = FALSE,
                                do.biexp = FALSE,
                                do.asinh = TRUE,
                                do.asym = TRUE,
                                asym_root = 2,
                                cofactor = 500,
                                stim_to_exclude = NULL)
{
  #- Checks
  if(is.null(stim_to_exclude) == FALSE){
    message(paste(stim_to_exclude, collapse = ", "), " stimulations are removed.")
  }

  #- Call extract_flow_exprs_data.R
  dt.exprs <- ICSR::extract_flow_exprs_data(gs = gs,
                                            output_nodes = output_nodes,
                                            parent_node = parent_node,
                                            cytokine_nodes = cytokine_nodes,
                                            do.comp = do.comp,
                                            do.biexp = do.biexp,
                                            do.asinh = do.asinh,
                                            do.asym = do.asym,
                                            asym_root = asym_root,
                                            cofactor = cofactor,
                                            stim_to_exclude = stim_to_exclude)

  #- Call extract_CYTNUM_data.R
  dt.cytnum <- ICSR::extract_CYTNUM_data(gs = gs,
                                         parent_node = parent_node,
                                         cytokine_nodes = cytokine_nodes,
                                         stim_to_exclude = stim_to_exclude)

  #- Output
  list("exprs" = dt.exprs,
       "cytnum" = dt.cytnum) %>%
    return()
}

