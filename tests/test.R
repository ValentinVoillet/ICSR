###--- R script
# CD4+ T
library(flowWorkspace)
library(tidyverse)
library(doMC)
library(here)
library(data.table)


# ml fhR/4.1.0-foss-2020b


#######################################
##############--- M72 ---##############
#######################################


####################################################
####################################################
###--- Open flowWorkspace objects
# /fh/fast/mcelrath_j/vvoillet_working/CHIL/Projects/CoVPN-3008/ICS
folders <- list.files(path = here("data-raw", "tmpdata"), full.names = TRUE)
gs_list <- foreach(i = folders) %do%
{
  cat(i, "\n")
  gs <- load_gs(path = i)
  gs %>% return()
}
dt <- enframe(gs_list, value = "flowWorkspace")
dt$name <- folders
dt <- dt %>%
  mutate(phenoData = map(.x = flowWorkspace, .f = function(x) return(x %>% pData)))


####################################################
####################################################
###--- Filtering - Remove duplicate samples
#- Duplicates
# Find duplicate samples
batchList <- lapply(1:nrow(dt), function(i)
{
  unique(pData(dt$flowWorkspace[[i]])[, c("BATCH", "PTID", "VISITNO", "STIM", "Run Num", "Replicate")])
})
batchDF <- do.call(rbind, batchList)
batchDF$`Run Num` <- as.numeric(as.character(batchDF$`Run Num`))
# Create ID tag for PTID:VISITNO:STIM
batchDF <- batchDF %>%
  mutate(PTID.VISITNO.STIM = paste(batchDF$PTID, batchDF$VISITNO, batchDF$STIM, sep = ":")) %>%
  arrange(PTID.VISITNO.STIM)
batchDF <- batchDF %>%
  group_by(PTID.VISITNO.STIM) %>%
  mutate(duplicated = ifelse(`Run Num` == max(`Run Num`), FALSE, TRUE))
batchDF$duplicated %>% table()

#- Filtering
bind_rows(dt$phenoData) %>%
  filter(!(PTID %in% "CTAC0100")) %>% # Remove control samples
  pull(sample_name) -> selected.samples
dt <- dt %>%
  mutate(filtering = map_chr(.x = phenoData, .f = function(x) ifelse(length(intersect(x$sample_name, selected.samples)) != 0, "Y", "N"))) %>%
  filter(filtering == "Y") %>%
  mutate(GatingSet.filtered = map2(.x = flowWorkspace, .y = phenoData, .f = function(x, y) x[which(y$sample_name %in% selected.samples), ] %>% return())) %>%
  mutate(phenoData.filtered = map(.x = GatingSet.filtered, .f = function(x) return(x %>% pData)))


####################################################
####################################################
###--- ICSR (package)
# usethis::create_github_token()
# credentials::set_github_pat("XXX") # put your token
# remotes::install_github("ValentinVoillet/ICSR", force = TRUE)


####################################################
####################################################
###--- extract_flow_exprs_data
parent_node <- "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+"
cytokine_nodes <- c("/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+/153+",
                    "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+/154+",
                    "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+/IFNg+",
                    "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+/IL17A_OR_IL17F",
                    "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+/IL2+",
                    "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+/IL4_OR_IL13",
                    "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+/GM-CSF+",
                    "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+/TNF+")
output_nodes <- c("/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+/R7+",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+/4+RA+",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+/Naive",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+/CM",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+/EM",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+/TEMRA",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+/153+",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+/154+",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+/IFNg+",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+/IL17A_OR_IL17F",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+/IL2+",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+/IL4_OR_IL13",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+/GM-CSF+",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+/TNF+",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+/CCR6+",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+/CXCR3+",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+/CXCR5+",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+/DR+",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+/Granulysin+",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+/Ki67+",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+/NKG2C+",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+/PD1+",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+/Perforin+")
lapply(X = dt$GatingSet.filtered, FUN = function(x) {
  ICSR::extract_flow_exprs_data(gs = x,
                                output_nodes = output_nodes,
                                parent_node = parent_node,
                                cytokine_nodes = cytokine_nodes,
                                do.comp = TRUE,
                                do.biexp = TRUE,
                                do.asinh = TRUE,
                                cofactor = 500) %>%
    return()}) %>%
  bind_rows() -> dt.results
dt.results_old <- readRDS(file = "/fh/fast/mcelrath_j/vvoillet_working/CHIL/Projects/M72_ID52-Frahm-Pilot/ICS/data/001_dt_CD4+T_1-of-8-cytokines(1).rds")


####################################################
####################################################
###--- extract_CYTNUM_data
parent_node <- "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+"
cytokine_nodes <- c("/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+/153+",
                    "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+/154+",
                    "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+/IFNg+",
                    "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+/IL17A_OR_IL17F",
                    "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+/IL2+",
                    "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+/IL4_OR_IL13",
                    "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+/GM-CSF+",
                    "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/S/L/3+/56-16-/gd-/4+/TNF+")
lapply(X = dt$GatingSet.filtered, FUN = function(x) {
  ICSR::extract_CYTNUM_data(gs = x,
                            parent_node = parent_node,
                            cytokine_nodes = cytokine_nodes) %>%
    return()}) %>%
  bind_rows() -> dt.results
dt.results_old <- readRDS(file = "/fh/fast/mcelrath_j/vvoillet_working/CHIL/Projects/M72_ID52-Frahm-Pilot/ICS/data/001_dt_CD4+T_1-of-8-cytokines(2).rds")


####################################################
####################################################
library(MIMOSA)

###--- runMIMOSA
inFileName <- here::here("output", "001_output", "in-mimosaSet.csv")
outFileName <- here::here("output", "001_output", "out-mimosaSet.csv")
ICSR::runMIMOSA(INFILE = inFileName, OUTFILE = outFileName)


####################################################
####################################################
###--- leiden_local
library(tidyverse)
library(here)
library(data.table)
library(doMC)
library(Rphenograph)
library(leiden)

#- Data
load(here("output", "001_output", "001_High-Dimensional-ICS-Analysis_CD4+T.RData"))

#- Leiden (asinh_asym)
markers <- colnames(dt)[str_detect(string = colnames(dt), pattern = "asinh_asym")]
partition.i <- ICSR::leiden_local(data = dt, markers = markers, k = 30, seed = 1234, niter = 100, res = 0.1)

