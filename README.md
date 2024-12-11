# singIST-reproducibility
![gabst_01dic](https://github.com/user-attachments/assets/0d36443d-007f-423b-ae75-64bb1a3e23c2)
This repository includes R scripts to reproduce all steps of singIST, including information from figures and tables from the original paper (bioRxiv, 2024). singIST is a method for comparative single-cell transcriptomics between disease models and humans, which provides explainable and quantitative insights into single-cell transcriptomics alignment. 

## Table of Contents

1. [Repository structure](#Repository-structure)
2. [Associated code of figures and tables](#Associated-code-of-figures-and-tables)
3. [Setup Instructions](#Setup-Instructions)
4. [Requirements](#Requirements)
5. [R Dependencies](#R-Dependencies)
   

# Repository structure
The repository is organized as follows:
- `0_rawdata`: human scRNA-seq [(Bangert, 2021)](https://pubmed.ncbi.nlm.nih.gov/33483337/) and disease models scRNA-seq Seurat objects.
- `1_input_preprocessing`: 
- 

Each report comes with a companion folder with the exported results. If the report name is 3_report.Rmd, the output folder will be 3_report_output/, so it is always caught by .gitignore. Important file locations are defined in the config.yml file and fecthed using the config R package.

# Associated code of figures and tables 
| Figure/Table      | R script                       |
| ------------ | --------------------------------- |
| Table 1 |  |
| Table 2  |               |
| Fig 3    | Texto en *cursiva*                |
| Fig 4    | Texto en *cursiva*                |
| Fig 5    | Texto en *cursiva*                |

# Setup Instructions
1. Clone the repository:
```bash
git clone https://github.com/amoruno/singIST-reproducibility.git
```
2. Download raw data:
```R
install.packages(c("Seurat", "googleCloudStorageR"))
library(googleCloudStorageR)
gcs_auth("path_to_your_service_account_key.json")
gcs_get_object("your_file.rds", bucket = "your_bucket_name", saveToDisk = "local_path/your_file.rds", overwrite = TRUE)
seurat_object <- readRDS("local_path/your_file.rds")

Sys.setenv("GCS_AUTH_FILE" = "path_to_your_service_account_key.json")
gcs_auth(Sys.getenv("GCS_AUTH_FILE"))
```
# Requirements
## System Requierements
### Databricks script
```
R version 4.2.2 (2022-10-31)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 22.04.5 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so

locale:
 [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
 [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
 [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
[10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
[1] SparkR_3.4.1   compiler_4.2.2 tools_4.2.2    Rserve_1.8-12 
```

### R script

# R Dependencies
