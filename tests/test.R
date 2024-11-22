###--- R script
library(flowCore)
library(flowWorkspace)
library(ICSR)
library(tidyverse)
library(here)
library(data.table)
library(doMC)


# ml R/4.4.0-gfbf-2023b



####################################################
####################################################
###--- ICSR (package)
# usethis::create_github_token()
# credentials::set_github_pat("XXX") # put your token
# remotes::install_github("ValentinVoillet/ICSR", force = TRUE)



###########################################################
##############--- CoVPN-3008 Case/Control ---##############
###########################################################


####################################################
####################################################
###--- Open flowWorkspace objects
# /fh/fast/mcelrath_j/vvoillet_working/CHIL/Projects/CoVPN-3008/ICS_Case-Control
folders <- list.files(path = here("data-raw", "tmpdata"), full.names = TRUE, pattern = "22")
gs_list <- foreach(i = folders) %do%
{
  cat(i, "\n")
  gs <- load_gs(path = i)
  pData(gs)$BATCH <- pData(gs)$Batch
  gs %>% return()
}


####################################################
####################################################
###--- compile_flow_events.R
parent_node <- "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+"
cytokine_nodes <- c("/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/IFNg+",
                    "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/IL2+",
                    "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/IL4+",
                    "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/IL5_OR_IL13+",
                    "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/IL17a+",
                    "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/IL21+",
                    "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/154+",
                    "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/TNFa+")
output_nodes <- c("/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/CCR7+",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/4+45RA+",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/N",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/CM",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/EM",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/TEMRA",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/IFNg+",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/IL2+",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/IL4+",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/IL5_OR_IL13+",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/IL17a+",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/IL21+",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/154+",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/TNFa+",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/GzB+",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/Perf+",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/25+",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/PD1+",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/CCR6+",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/CXCR3+",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/CXCR5+",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/Ki67+",
                  "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/FoxP3+")
#- Processing
list.res <- lapply(X = gs_list, FUN = function(x) {
  message(pData(x)$BATCH %>% unique())
  compile_flow_events(gs = x,
                      output_nodes = output_nodes,
                      parent_node = parent_node,
                      cytokine_nodes = cytokine_nodes,
                      do.comp = FALSE,
                      do.biexp = FALSE,
                      do.asinh = TRUE,
                      do.asym = TRUE,
                      cofactor = 500,
                      stim_to_exclude = "sebctrl") %>%
    return()})

#- Checks
# list.res[[1]]$exprs %>% head()
# list.res[[1]]$cytnum %>% head()

#- Output
saveRDS(object = list.res, file = here("data", "raw_test.rds"))


####################################################
####################################################
###--- Filtering
#- Unreliable
unreliable.s <- readxl::read_xlsx("~/Documents/CHIL/CoVPN-3008/ICS_Case-Control/misc-docs/CoVPN 3008 Case Cohort_Unreliable data_16NOV23_A_modified_08FEB2024.xlsx")
unreliable.s <- unreliable.s[2:nrow(unreliable.s), ]
unreliable.s %>%
  filter(`Unreliable Stim` == "ALL") %>%
  mutate(SAMPLE = paste(`Batch #`, PTID, VISITNO, `Run Num`, sep = "_")) %>%
  pull(SAMPLE) -> unreliable.s_1
unreliable.s %>%
  filter(!(`Unreliable Stim` %in% c("ALL", "negctrl"))) %>%
  mutate(SAMPLE = paste(`Batch #`, PTID, VISITNO, `Run Num`, "1", `Unreliable Stim`, sep = "_")) %>%
  pull(SAMPLE) -> unreliable.s_2

#- dt
list.res.raw <- readRDS(file = "tests/data/raw_test.rds")
lapply(X = list.res.raw, function(x) x$exprs) %>%
  bind_rows() -> dt.exprs.raw
lapply(X = list.res.raw, function(x) x$cytnum) %>%
  bind_rows() -> dt.cytnum.raw

#- Filtering - dt.exprs
dt.MTL <- readxl::read_xlsx("~/Documents/CHIL/CoVPN-3008/ICS_Case-Control/misc-docs/CoVPN3008 Cases_MTL_A_14Sep23_VV.xlsx")
colnames(dt.MTL)[1] <- "BATCH"
colnames(dt.MTL)[10] <- "RUNNUM"
do_filtering(dt = dt.exprs.raw,
             remove_dup = TRUE,
             dt.MTL_dup = dt.MTL,
             remove_participants = TRUE,
             participants_list = "CTAC0036",
             remove_nsub = FALSE,
             cutoff_nsub = NULL,
             flag_unreliable = TRUE,
             unreliable_list_all = unreliable.s_1,
             unreliable_list_by_stim = unreliable.s_2) -> dt.exprs
dt.exprs$reliable_flag %>% table() # 504 & 195069
dim(dt_old_filt) # 195069
dt.exprs %>%
  filter(reliable_flag == 1) -> dt.exprs
identical(paste(dt.exprs$FCS, dt.exprs$`/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/25+`),
          paste(dt_old_filt$FCS, dt_old_filt$`/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/25+`))  # TRUE

#- Filtering - dt.cytnum
do_filtering(dt = dt.cytnum.raw,
             remove_dup = TRUE,
             dt.MTL_dup = dt.MTL,
             remove_participants = TRUE,
             participants_list = "CTAC0036",
             remove_nsub = FALSE,
             cutoff_nsub = NULL,
             flag_unreliable = TRUE,
             unreliable_list_all = unreliable.s_1,
             unreliable_list_by_stim = unreliable.s_2) -> dt.cytnum


####################################################
####################################################
###--- QC
create_report_QC_ICS(dt.exprs = dt.exprs,
                     dt.cytnum = dt.cytnum,
                     output_dir = "tests/",
                     report_title = "CHIL ICS QC Report",
                     report_author = "Valentin Voillet")



