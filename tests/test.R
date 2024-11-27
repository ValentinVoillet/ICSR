###--- R script
library(flowCore)
library(flowWorkspace)
#library(ICSR)
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
ICSR::do_filtering(dt = dt.exprs.raw,
             remove_dup = TRUE,
             dt.MTL_dup = dt.MTL,
             remove_participants = TRUE,
             participants_list = "CTAC0036",
             remove_nsub = FALSE,
             cutoff_nsub = NULL,
             flag_unreliable = TRUE,
             unreliable_list_all = unreliable.s_1,
             unreliable_list_by_stim = unreliable.s_2) -> dt.exprs
dt.exprs$reliable_flag %>% table() # 504 & 195069 # Several check were made.

#- Filtering - dt.cytnum
ICSR::do_filtering(dt = dt.cytnum.raw,
             remove_dup = TRUE,
             dt.MTL_dup = dt.MTL,
             remove_participants = TRUE,
             participants_list = "CTAC0036",
             remove_nsub = FALSE,
             cutoff_nsub = NULL,
             flag_unreliable = TRUE,
             unreliable_list_all = unreliable.s_1,
             unreliable_list_by_stim = unreliable.s_2) -> dt.cytnum

#- Output
dt.exprs %>% saveRDS(file = "tests/data/dt_exprs.rds")
dt.cytnum %>% saveRDS(file = "tests/data/dt_cytnum.rds")


####################################################
####################################################
###--- QC
dt.exprs <- readRDS(file = "tests/data/dt_exprs.rds")
dt.cytnum <- readRDS(file = "tests/data/dt_cytnum.rds")

cytokine_nodes <- c("/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/IFNg+",
                    "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/IL2+",
                    "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/IL4+",
                    "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/IL5_OR_IL13+",
                    "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/IL17a+",
                    "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/IL21+",
                    "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/154+",
                    "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14-/L/S/19-/3+/3+excl 16br/56-16-/4+/TNFa+")
markers <- c("asinh_asym_CD4", "asinh_CD8", "asinh_asym_CD8")
create_report_QC_ICS(dt.exprs = dt.exprs,
                     dt.cytnum = dt.cytnum,
                     cytokine_nodes = cytokine_nodes,
                     markers = markers,
                     output_dir = "tests/",
                     report_title = "CHIL ICS QC Report",
                     report_author = "Valentin Voillet")


####################################################
####################################################
###--- Leiden
dt.test <- dt.exprs %>%
  filter(reliable_flag == 1)
set.seed(1234)
n <- sample(1:nrow(dt.test), 10000)
markers <- colnames(dt.test)[str_detect(string = colnames(dt.test), pattern = "asinh_asym")]
# use_condaenv("/opt/miniconda3/envs/spyder-env/bin/python", required = TRUE)
x = dt.test[n, ]
leiden.partitionx <- ICSR::leiden_local(data = x, markers = markers[1:6], k = 30, res = 1, niter = 2, seed = 1234)
leiden.partition <- leiden_local(data = x, markers = markers[1:6], k = 30, res = 1, niter = 2, seed = 1234)
identical(leiden.partition, leiden.partitionx)
x$LEIDEN <- factor(paste("Leiden:", leiden.partition))


####################################################
####################################################
###--- MIMOSA
#- Process as usual
x$Leiden <- factor(paste("Leiden:", leiden.partition))
dt.summary <- x %>%
  group_by(BATCH, PTID, STIM, VISITNO, RUNNUM, NSUB) %>%
  summarize(CYTNUM = n()) %>%
  group_by(BATCH, PTID, STIM, VISITNO, RUNNUM) %>%
  summarise(NSUB = sum(NSUB), CYTNUM = sum(CYTNUM)) %>%
  ungroup()

#- Pre-processing by Leiden
dt.tmp <- x %>%
  group_by(BATCH, PTID, STIM, VISITNO, Leiden, .drop = FALSE) %>%
  summarize(CYTNUM = n())
dt.tmp <- dt.tmp %>%
  mutate(NSUB = plyr::mapvalues(x = paste(BATCH, PTID, STIM, VISITNO),
                                from = paste(dt.summary$BATCH,
                                             dt.summary$PTID,
                                             dt.summary$STIM,
                                             dt.summary$VISITNO),
                                to = dt.summary$NSUB,
                                warn_missing = FALSE)) %>%
  mutate(NSUB = NSUB %>% as.numeric()) %>%
  ungroup() %>%
  mutate(SAMPLE = paste(PTID, VISITNO, Leiden)) %>%
  select(BATCH, PTID, STIM, VISITNO, SAMPLE, Leiden, NSUB, CYTNUM)
#- Background subtraction
table(dt.tmp$STIM)
dt.tmp_1 <- dt.tmp %>%
  filter(STIM == "negctrl") %>%
  rename(NSUB_NEG = "NSUB", CYTNUM_NEG = "CYTNUM") %>%
  select(SAMPLE, NSUB_NEG, CYTNUM_NEG)
dt.tmp_2 <- dt.tmp %>%
  filter(STIM != "negctrl")
dt.Leiden <- merge(x = dt.tmp_2, y = dt.tmp_1, by = "SAMPLE", all.x = TRUE) %>%
  mutate(PCTPOS = (CYTNUM / NSUB) * 100) %>%
  mutate(PCTNEG = (CYTNUM_NEG / NSUB_NEG) * 100) %>%
  mutate(PCTPOS_ADJ = PCTPOS - PCTNEG) %>%
  select(BATCH, PTID, STIM, VISITNO, Leiden, NSUB, CYTNUM, PCTPOS, NSUB_NEG, CYTNUM_NEG, PCTNEG, PCTPOS_ADJ) %>%
  arrange(BATCH, PTID, STIM, VISITNO)
dt.Leiden <- na.omit(dt.Leiden)

#- Data
mimosaSet <- dt.Leiden %>%
  mutate(SAMPLE = paste(PTID, VISITNO)) %>%
  select(PTID, STIM, VISITNO, SAMPLE, Leiden, NSUB, CYTNUM, NSUB_NEG, CYTNUM_NEG)
names(mimosaSet) <- toupper(names(mimosaSet))
write.table(x = mimosaSet, file = here("tests", "data", "in-mimosaSet.csv"), row.names = FALSE, sep = ",")

#- Run MIMOSA
# A multiplicity adjustment is made to the individual peptide pool p-values
inFileName <- here::here("tests", "data", "in-mimosaSet.csv")
outFileName <- here::here("tests", "data", "out-mimosaSet.csv")
ICSR::runMIMOSA(INFILE = inFileName, OUTFILE = outFileName, CLUSTERS = paste("Leiden:", 1))


####################################################
####################################################
###--- UMAP
umap_test <- UMAP_local(dt = dt.test, markers = markers)


####################################################
####################################################
###--- annotation
dt <- dt.test
dt$Leiden <- sample(paste("Leiden:", 1:8), nrow(dt), replace = T)

# colnames
colnames(dt)
colnames(dt)[11:34] <- str_remove(string = colnames(dt)[11:34], pattern = "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14[-]/L/S/19[-]/3[+]/3[+]excl 16br/56[-]16[-]/4[+]/")
colnames(dt)[13] <- "CD45RA+"
colnames(dt)[11] <- "CD4+"
colnames(dt)[28] <- "CD25+"
colnames(dt)[24] <- "CD154+"
colnames(dt)[25] <- "TNF+"
colnames(dt)[14:17] <- paste0(colnames(dt)[14:17], "+")

flowjo_gates <- colnames(dt)[str_detect(string = colnames(dt), pattern = "[+]")]
flowjo_gates <- flowjo_gates[-c(22)]

colnames(dt)[62] <- "asinh_asym_Perf"
colnames(dt)[63] <- "asinh_asym_IL5_OR_IL13"
colnames(dt)[67] <- "asinh_asym_GzB"
colnames(dt)[74] <- "asinh_asym_PD1"
colnames(dt)[87] <- "asinh_asym_FoxP3"
markers <- colnames(dt)[str_detect(string = colnames(dt), pattern = "asinh_asym")]
markers <- markers[-c(9, 10, 14, 17, 19, 20, 22, 23)]
markers

cluster_col <- "Leiden"

# dotplot
order_level_marker <- c("N", "CM", "EM", "TEMRA", "CD4", str_remove_all(string = markers, pattern = "asinh_asym_"))
dotplot_annotation(dt = dt,
                   markers = markers,
                   flowjo_gates = flowjo_gates,
                   cluster_col = "Leiden",
                   level_marker = order_level_marker,
                   order = paste("Leiden:", c(8, 6, 7, 1:5)))

# M72
load("~/Documents/CHIL/M72_ID52-Frahm-Pilot/ICS/output/001_output/001_High-Dimensional-ICS-Analysis_CD4+T.RData")
colnames(dt)
colnames(dt)[15:38] <- str_remove(string = colnames(dt)[15:38], pattern = "/Time/K1/K2/K3/K4/K5/K6/K7/K8/Lv/14[-]/S/L/3[+]/56[-]16[-]/gd[-]/4[+]/")
colnames(dt)[15] <- "CD4+"
colnames(dt)[16] <- "CCR7+"
colnames(dt)[17] <- "CD45RA+"
colnames(dt)[22] <- "CD153+"
colnames(dt)[23] <- "CD154+"
colnames(dt)[25] <- "IL17a_IL17F+"
colnames(dt)[27] <- "IL4_IL13+"
colnames(dt)[28] <- "GMCSF+"
colnames(dt)[33] <- "HLADR+"
colnames(dt)[34] <- "GranLys+"
colnames(dt)[18:21] <- paste0(colnames(dt)[18:21], "+")
flowjo_gates <- colnames(dt)[str_detect(string = colnames(dt), pattern = "[+]")]
markers <- colnames(dt)[str_detect(string = colnames(dt), pattern = "asinh_asym")]
order_level_marker <- c("Naive", "CM", "EM", "TEMRA", "CD4", str_remove_all(string = markers, pattern = "asinh_asym_"))
ICSR::dotplot_annotation(dt = dt,
                   markers = markers,
                   flowjo_gates = flowjo_gates,
                   cluster_col = "Leiden",
                   level_marker = order_level_marker,
                   order = paste("Leiden:", c(7, 10, 5, 9, 1, 2, 6, 20, 12, 19, 13, 14, 3, 4, 8, 15, 11, 16, 18, 17, 21)))

