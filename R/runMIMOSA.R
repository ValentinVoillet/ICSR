#' Run MIMOSA Response Calling
#'
#' This function applies the MIMOSA (Mixture Models for Single-Cell Analysis) framework
#' to identify significant vaccine-induced responses. It models the distribution of
#' cytokine-positive cells in stimulated vs. unstimulated samples to determine
#' "responder" status per subject and per cluster.
#'
#' @param INFILE Character. Path to the input .csv file containing aggregated counts
#'   (must include \code{NSUB}, \code{CYTNUM}, \code{PTID}, \code{STIM}, \code{VISITNO}, and \code{LEIDEN}).
#' @param OUTFILE Character. Path where the resulting MIMOSA calls will be saved as a .csv.
#' @param CLUSTERS Character vector. The cluster labels to iterate over (e.g., \code{paste0("Cluster_", 1:10)}).
#' @param MIMOSA_THRESHOLD_FDR Numeric. False Discovery Rate (FDR) threshold for response calling. Default is 0.01.
#' @param FIT_METHOD Character. Method for model fitting: \code{"mcmc"} (Markov Chain Monte Carlo) or \code{"EM"} (Expectation-Maximization).
#' @param COMBINE Character vector. Variables to combine during stratification (e.g., \code{"VISITNO"}).
#' @param ANTIGENFILTER Character vector. List of antigens (STIM) to exclude from the analysis.
#' @param MINSAMPLES Integer. Minimum number of samples required. Default is 10.
#' @param ITER Integer. Number of MCMC iterations. Minimum is 250,000.
#' @param BURN Integer. Number of burn-in iterations for MCMC. Minimum is 50,000.
#'
#' @return A \code{.csv} file saved to \code{OUTFILE} containing probabilities of response and FDR-adjusted calls.
#'
#' @import MIMOSA tidyverse data.table plyr
#' @export
runMIMOSA <- function(INFILE = NULL,
                      OUTFILE = NULL,
                      CLUSTERS = paste("Leiden:", 1:10),
                      MIMOSA_THRESHOLD_FDR = 0.01,
                      FIT_METHOD = "mcmc",
                      COMBINE = NULL,
                      ANTIGENFILTER = NULL,
                      MINSAMPLES = 10,
                      ITER = 250000,
                      BURN = 50000)
{
  # 1. Argument Validation & Housekeeping
  # ---------------------------------------------------------------------------
  if(!(exists("INFILE") & exists("OUTFILE"))) { stop("<infile.csv> or <oufile.csv> not defined\n") }
  if(exists("FIT_METHOD")) { match.arg(FIT_METHOD,c("EM", "mcmc")) }
  if(exists("COMBINE")){ COMBINE <- COMBINE[COMBINE %in% c("STIM", "VISITNO")] }
  message("Will filter for antigens: ", ANTIGENFILTER)
  message("fitmethod is ", FIT_METHOD)
  if(!is.null(COMBINE)) message("combining observations for ", COMBINE)
  message("Adjustment parameter is currently ignored. Adjusting across antigens by default.")
  if(!(as.numeric(ITER) > 250000)){
    ITER <- 250000
  } else {
    ITER <- as.numeric(ITER)
  }
  if(!(as.numeric(BURN) > 50000)){
    BURN <- 50000
  } else {
    BURN <- as.numeric(BURN)
  }

  # 2. Data Preparation
  # ---------------------------------------------------------------------------
  # Open
  data <- read.csv(INFILE)
  colnames(data) <- toupper(colnames(data))
  if(!is.null(ANTIGENFILTER)){
    message("Applying antigen filter ", ANTIGENFILTER)
    data <- subset(data, !ANTIGEN %in% ANTIGENFILTER)
  }
  # Ensure critical variables are factors for MIMOSA ExpressionSet construction
  tofactors <- intersect(c("PTID", "STIM", "VISITNO", "SAMPLE", "LEIDEN"), colnames(data))
  data[tofactors] <- lapply(data[tofactors], factor)
  # Check if the results columns are present
  if(!all(c("NSUB", "CYTNUM") %in% colnames(data))){
    stop("NSUB or CYTNUM columns are not present. Stopping")
    q(save = "no",status = 1)
  }

  # 3. Construct MIMOSA ExpressionSet
  # ---------------------------------------------------------------------------
  # MIMOSA requires a specific ExpressionSet object where 'reference' is the unstimulated control
  annotations <- c("PTID", "STIM", "VISITNO", "LEIDEN")
  F <- as.formula(paste("component ~", paste(annotations, collapse = "+"), sep = ""))
  # Construct the data
  E <- ConstructMIMOSAExpressionSet(thisdata = data,
                                    reference = NULL,
                                    measure.columns = c("CYTNUM", "NSUB"),
                                    other.annotations = annotations,
                                    default.cast.formula = F,
                                    .variables = .c(PTID, STIM, VISITNO), # The .variables argument describes the grouping variables used to stratify the data into a single experimental unit. Within each unique set of levels in that group, we'll have all the different leiden clusters for a subject.
                                    featureCols = 1,
                                    ref.append.replace = "_NEG")

  # 4. Stratification Logic
  # ---------------------------------------------------------------------------
  # Dynamically build the formula based on stratification variables
  groupBy <- function(var, F)
  {
    F <- gsub(paste("\\+", var, "\\+", sep=""), "\\+", gsub(" ", "", paste(deparse(F), collapse = "")))
    if(grepl("\\|", F)){
      return(as.formula(gsub("(\\|.*)", paste("\\1+", var, sep=""), F)))
    }else{
      return(as.formula(gsub("$", paste("|", var, sep=""), F)))
    }
  }
  componentY <- function(F)
  {
    F <- gsub(" ", "", paste(deparse(F), collapse = ""))
    as.formula(gsub("component(*~)", "NSUB+CYTNUM\\1",F))
  }
  STRATIFY <- c("VISITNO", "STIM")
  STRATIFY <- setdiff(STRATIFY, COMBINE)
  for(i in STRATIFY){
    F <- groupBy(i, F)
  }
  F <- componentY(F)
  F <- list(F)
  # F[[1]] <- as.formula("NSUB + CYTNUM ~ PTID + LEIDEN | STIM + VISITNO") # == NSUB + CYTNUM ~ PTID | STIM + VISITNO
  # In this case we fit a separate model for each visit and stimulation across all subjects.
  # Below we loop on leiden, so that's also a separate model per leiden.

  # 5. Model Fitting (Loop over Clusters)
  # ---------------------------------------------------------------------------
  # NOTE: called separately for each cluster (Daryl: it was done for each cytokine (IL2_or_IFNg, etc..)).
  # Here, we consider that clusters are different combinations of cytokines/markers.
  # Each element of result is a MIMOSAResultList (a list of MIMOSAResult's) the MIMOSAResultList has one element for each combo of stratification vars
  result <- NULL
  for(ck in CLUSTERS)
  {
    Esub <- E[, E$LEIDEN == ck]
    result <- c(result,
                lapply(F, function(F) MIMOSA(F,
                                             Esub,
                                             ref = RefTreat %in% "Reference",
                                             subset = RefTreat %in% "Treatment",
                                             method = FIT_METHOD,
                                             iter = ITER,
                                             burn = BURN)))
  }

  # 6. Result Aggregation & FDR Correction
  # ---------------------------------------------------------------------------
  # Flatten the nested result lists into a single data frame
  res.full <- ddply(do.call(rbind, lapply(result, function(r) do.call(rbind, lapply(r, function(tc) {
    cbind(pData(tc),
          effect = {
            if(class(tc@result) == "MDMixResult"){
              prop.table(as.matrix(tc@result@data$n.stim), 1)[,2] -
                prop.table(as.matrix(tc@result@data$n.unstim), 1)[,2]
            }else{
              prop.table(as.matrix(tc@result@n.stim),1)[,2] -
                prop.table(as.matrix(tc@result@n.unstim),1)[,2]
            }
          },
          Pr_nonresp = tc@z[,1],
          Pr_resp = tc@z[,2]
    )})))),
    .(STIM), transform, # == group_by(STIM)
    fdr = MIMOSA:::fdr(cbind(Pr_nonresp, Pr_resp)))
  res.full$Mimosa.call <- res.full$fdr < MIMOSA_THRESHOLD_FDR

  # 7. Output
  # ---------------------------------------------------------------------------
  with(res.full, if(all(Pr_nonresp == fdr)){
    warning("The FDR is equal to the probability of non-response.\n",
            "Unless you are adjusting over one antigen something may have gone wrong.");
  })
  colnames(res.full) <- tolower(colnames(res.full))
  write.table(res.full, file = OUTFILE, sep = ",", row.names = FALSE, quote = FALSE)
}

