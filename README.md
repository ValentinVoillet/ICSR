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
