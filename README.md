# ICSR

## 0. Getting Started

To start, you need:

-   **Master thaw list**: A file containing both `Batch` and `SAMP_ORD` columns to define your experimental structure and processing order.
-   **.CSV file**: A `.csv` file mapping your raw data files (`.fcs` and `.xml`) to specific batches (assayid).

**Manifest Example**

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
