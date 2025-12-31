# ICSR

## 0. Getting Started

To start, you need:

-   **Master thaw list**: A file containing both `Batch` and `SAMP_ORD` columns to define your experimental structure and processing order.
-   **.CSV file**: A `.csv` file mapping your raw data files (`.fcs` and `.xml`) to specific batches (assayid).

**Example of .CSV file**

| assayid | xml                  | fcs                    |
|:--------|:---------------------|:-----------------------|
| ID_882  | paths/to/file_01.xml | paths/to/sample_01.fcs |
| ID_883  | paths/to/file_01.xml | paths/to/sample_02.fcs |
| ID_884  | paths/to/file_02.xml | paths/to/sample_03.fcs |

------------------------------------------------------------------------

## Step 1: Creation of GatingSet Object

The first step in the `ICSR` workflow is the creation of a **GatingSet** object.

A **GatingSet** is a powerful data structure provided by the `flowWorkspace` ecosystem. Think of it as a specialized container that holds not just your raw single-cell data (FCS files), but also the entire "gating hierarchy" (the XML Workspace). By using a GatingSet, `ICSR` can efficiently manage thousands of gates across multiple batches while ensuring statistical consistency.

**Important:** The `GatingSet` objects are saved locally in a folder named `data-raw/tmpdata/`. You must create this folder directory before running the extraction function.

You can initialize your GatingSet using the following code:

```{r}
library(flowCore)
library(flowWorkspace)
library(CytoML)
library(tidyverse)
library(here)
library(data.table)
library(parallel)
library(ICSR)

#- Open Master Thaw List & .CSV file
batchData <- readxl::read_xlsx(path = "MTL.xlsx")
dt.workspace <- read_csv(file = "batch_to_workspace_map.csv")

#- GatingSet Obj.
cl <- makeCluster(5)
clusterExport(cl, c("batchData", "dt.workspace", "create_gating_set"))
parApply(cl, dt.workspace, 1, function(x) {
  ICSR::raw_to_GatingSet(
    assayid = x[["assayid"]],
    xml_path = x[["xml"]],
    fcs_path = x[["fcs"]]
  )
})
stopCluster(cl)
```

------------------------------------------------------------------------

## Step 2: Extract Data of Interest

Once your **GatingSet** objects are created, you can extract specific biological data for downstream analysis. The function `compile_flow_events` allows you to pull two types of information simultaneously:

1.  **Fluorescence Intensities:** The actual signal strength for every marker on a per-cell basis (supporting both raw and transformed values).

**Example of `head(output$exprs)`:**

| BATCH | PTID      | STIM | RUNNUM | REPLICATE | CYTNUM | asinh_IFNg | asinh_IL2 | IFNg+ | IL2+  |
|:-------|:-------|:-------|:-------|:-------|:-------|:-------|:-------|:-------|:-------|
| 2410  | Subject_A | Env  | 1      | 1         | 2      | 4.52       | 3.12      | TRUE  | TRUE  |
| 2410  | Subject_A | Env  | 1      | 1         | 1      | 0.82       | 4.05      | FALSE | TRUE  |
| 2410  | Subject_B | Env  | 2      | 1         | 0      | 0.15       | 0.45      | FALSE | FALSE |

2.  **Cell Counts:** Population statistics for specific gates across your entire experiment.

**Example of `head(output$cytnum)`:**

| BATCH | PTID      | STIM | RUNNUM | REPLICATE | NSUB  | boolean_CYTNUM | CYTNUM |
|:------|:----------|:-----|:-------|:----------|:------|:---------------|:-------|
| 2410  | Subject_A | Env  | 1      | 1         | 45000 | TRUE           | 152    |
| 2410  | Subject_A | Env  | 1      | 1         | 45000 | FALSE          | 44848  |
| 2410  | Subject_B | Env  | 2      | 1         | 42300 | TRUE           | 45     |
| 2410  | Subject_B | Env  | 2      | 1         | 42300 | FALSE          | 42255  |

### Example: Extracting Cytokine-Positive CD4+ T Cells

In the following code, we focus on a "Boolean" population—specifically, CD4+ T cells expressing at least one cytokine (CD153, CD154, IFNg, IL17A/F, IL2, IL4/13, GM-CSF and/or TNF).

```{r}
#- Open GatingSet obj.
folders <- list.files(path = here("data-raw", "tmpdata"), full.names = TRUE)
gs_list <- lapply(folders, load_gs)

#- Variables
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

#- Extraction
list.res <- lapply(X = gs_list, FUN = function(x) {
  message(pData(x)$BATCH %>% unique())
  compile_flow_events(
    gs = x,
    output_nodes = output_nodes,
    parent_node = parent_node,
    cytokine_nodes = cytokine_nodes,
    pData_cols = c("BATCH", "SAMP_ORD", "PTID", "STIM", "VISITNO", "Run Num", "Collection Num", "Replicate"),
    do.comp = FALSE,
    do.biexp = FALSE,
    do.asinh = TRUE,
    do.asym = TRUE,
    cofactor = 500,
    stim_to_exclude = "sebctrl"
    ) %>%
   return()
  })
```

------------------------------------------------------------------------

## Step 3: Data Filtering & Quality Control

Before proceeding to clustering or statistical inference, it is essential to clean your dataset. Analyzing samples with insufficient cell counts or unreliable data can lead to false positives and batch effect artifacts.

In `ICSR`, we recommend a standard filtering workflow:

1.  **Remove Control/Unreliable Samples:** Drop technical controls or samples flagged during wet-lab processing.
2.  **Handle Duplicates:** Ensure only the most reliable replicate is kept if necessary.
3.  **Minimum Event Threshold:** Depending on your study, filter out samples where the parent population count (`NSUB`) is too low to provide statistical power (e.g., needing at least 10,000 CD4+ T cells).

### Example Filtering Workflow

```{r}
#- Filtering
dt.exprs <- list.res$exprs %>%
  bind_rows() %>%
  filter(!(PTID %in% "CTACX"))
dt.cytnum <- list.res$cytnum %>%
  bind_rows() %>%
  filter(!(PTID %in% "CTACX"))
```

------------------------------------------------------------------------

## Step 4: Quality Control Reporting

After filtering and transformation, it is critical to perform a final inspection of the data. The `create_report_QC_ICS` function generates an automated, interactive HTML report that summarizes the state of your experiment.

This report allows you to: \* **Validate Transformations:** Confirm that the Arcsinh cofactors successfully resolved "squashing" or negative value artifacts. \* **Inspect Marker Distributions:** Use ridge plots to check for consistent staining across different batches. \* **Review Polyfunctionality:** Visualize the distribution of `CYTNUM` (number of cytokines per cell) across stimulations.

### Example Clustering Workflow

``` r
library(rmarkdown)

#- Define the markers you want to inspect in the ridge plots
my_markers <- c("asinh_IFNg", "asinh_IL2", "asinh_TNF", "asinh_CD154", "asinh_CD107a")

#- Create the automated HTML report
create_report_QC_ICS(
  dt.exprs = dt.exprs,
  dt.cytnum = dt.cytnum,
  markers = my_markers,
  cytokine_nodes = c("IFNg", "IL2", "TNF"),
  report_title = "Post-Filtering QC",
  report_author = "Your Name / Lab"
)
```

------------------------------------------------------------------------

## Step 5: High-Dimensional Clustering

Once the data is filtered and transformed, we identify cell populations using unsupervised clustering. Different clustering methods can be used. Here, we run the `leiden_local` function which applies the Leiden algorithm to the single-cell data.

1.  **Select Markers:** Choose phenotypic markers to define subsets (e.g., Memory vs Naive).
2.  **Run Leiden:** Execute the algorithm to partition the cells into clusters.
3.  **Assign Clusters:** Add the cluster labels back to your main dataset for downstream analysis.

### Example Clustering Workflow

```{r}
#- Define markers for clustering
markers <- colnames(dt.exprs)[str_detect(string = colnames(dt.exprs), pattern = "asinh_asym")]

#- Run Leiden clustering
partition <- leiden_local(data = dt.exprs, markers = markers, k = 30, res = 1, niter = 1, seed = 1234)
dt.exprs$LEIDEN <- paste("LEIDEN", partition, sep = "_")

#- Preview results
head(dt.exprs[, .(PTID, STIM, cluster)])
```

------------------------------------------------------------------------

## Step 6: UMAP Dimensionality Reduction

To visualize the high-dimensional data and the identified clusters in a 2D space, we use the **UMAP** (Uniform Manifold Approximation and Projection) algorithm. This allows for a visual validation of the clustering results and the identification of spatial relationships between populations.

1.  **Subsampling:** Since UMAP is computationally intensive, it is common to run it on a representative subset of cells (e.g., 50,000 cells).
2.  **Select Markers:** Use the same phenotypic markers used for clustering to ensure the 2D map reflects the same biological backbone.
3.  **Execution:** Run the UMAP algorithm and join the coordinates back to your data for plotting.

### Example UMAP Workflow

```{r}
#- Define markers for clustering
markers <- colnames(dt.exprs)[str_detect(string = colnames(dt.exprs), pattern = "asinh_asym")]

#- Run UMAP
UMAP.emb <- ICSR::UMAP_local(dt = dt.exprs, markers = markers, n_neighbors = 10, min_dist = .1, verbose = TRUE, n_components = 2, seed = 1234)
dt.exprs$UMAP_1 <- UMAP.emb[, 1]
dt.exprs$UMAP_2 <- UMAP.emb[, 2]

#- Preview results
head(dt.exprs[, .(PTID, STIM, UMAP_1, UMAP_2, cluster)])
```

------------------------------------------------------------------------

## Step 7: Statistical Response Calling with MIMOSA

To determine if the observed cytokine production in a cluster is a true biological response or just background noise, we use the `runMIMOSA` function. This applies a Bayesian framework to compare stimulated samples against their respective negative controls.

1.  **Aggregated Input:** MIMOSA requires a table of counts (e.g., the output of your cluster-level summary) containing `NSUB` (total cells) and `CYTNUM` (positive cells).
2.  **Cluster Iteration:** The function fits a model for each identified phenotype (Leiden clusters) to see which specific subsets are responding.
3.  **FDR Correction:** It automatically calculates the False Discovery Rate (FDR) across antigens to provide a robust "Response Call" (TRUE/FALSE).

### Example MIMOSA Workflow

This step is quite code-heavy because of the data wrangling required for MIMOSA (joining negative controls).

```{r}
#- Summary
dt.summary <- dt.exprs %>%
  group_by(BATCH, PTID, STIM, VISITNO, RUNNUM, NSUB) %>%
  summarize(CYTNUM = n()) %>%
  group_by(BATCH, PTID, STIM, VISITNO, RUNNUM) %>%
  summarise(NSUB = sum(NSUB), CYTNUM = sum(CYTNUM)) %>%
  ungroup()

#- Pre-processing by LEIDEN
dt.tmp <- dt.exprs %>%
  group_by(BATCH, PTID, STIM, VISITNO, LEIDEN, .drop = FALSE) %>%
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
  mutate(SAMPLE = paste(PTID, VISITNO, LEIDEN)) %>%
  select(BATCH, PTID, STIM, VISITNO, SAMPLE, LEIDEN, NSUB, CYTNUM)

#- Background subtraction
table(dt.tmp$STIM) # negctrl or NEGCTRL
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
  select(BATCH, PTID, STIM, VISITNO, LEIDEN, NSUB, CYTNUM, PCTPOS, NSUB_NEG, CYTNUM_NEG, PCTNEG, PCTPOS_ADJ) %>%
  arrange(BATCH, PTID, STIM, VISITNO)
dt.Leiden <- na.omit(dt.Leiden)

#- .CSV
mimosaSet <- dt.Leiden %>%
  mutate(SAMPLE = paste(PTID, VISITNO)) %>%
  select(PTID, STIM, VISITNO, SAMPLE, LEIDEN, NSUB, CYTNUM, NSUB_NEG, CYTNUM_NEG)
names(mimosaSet) <- toupper(names(mimosaSet))
write.table(x = mimosaSet, file = "MIMOSA_in.csv", row.names = FALSE, sep = ",")

#- Run MIMOSA response calling
# This will save a CSV with probabilities (Pr_resp) and FDR calls
runMIMOSA(
  INFILE = "MIMOSA_in.csv",
  OUTFILE = "MIMOSA_out.csv",
  CLUSTERS = paste0("LEIDEN_", 1:10),
  MIMOSA_THRESHOLD_FDR = 0.01,
  FIT_METHOD = "mcmc"
)

#- Load results to see responders
mimosa_res <- read.csv("MIMOSA_out.csv")
head(mimosa_res[mimosa_res$mimosa_call == TRUE, ])
```

------------------------------------------------------------------------

## Step 8: Downstream Analysis & Visualization

With the clusters identified and the statistical response calls completed, you can now perform deep-dive analyses into your specific biological questions. Since the `dt.exprs` and `mimosa_res` objects contain all the necessary metadata, the visualization possibilities are vast and depend on your study's objectives.

Common downstream tasks include:

1.  **Phenotype Mapping:** Using boxplots or heatmaps to visualize the expression of markers across clusters to "name" the populations (e.g., identifying Cluster 5 as Th1-like).
2.  **Response Landscape:** Visualizing the percentage of responding cells per cluster across different clinical groups or timepoints.
3.  **Spatial Distribution:** Projecting clusters back onto the UMAP to see where the clusters are physically located in the high-dimensional space.
4.  **Polyfunctionality:** Analyzing the combinations of cytokines produced within each specific cluster.

### Example Visualization Concepts

While the analysis is user-defined, common outputs at this stage include:

-   **UMAP Plots:** Colored by Cluster, PTID, or BATCH to identify patterns.
-   **Frequency Plots:** Boxplots showing the percentage of a specific cluster relative to the parent population.
-   **MFI Heatmaps:** Summarizing marker expression per cluster for easy phenotype identification.

