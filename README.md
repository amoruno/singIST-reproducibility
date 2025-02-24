# singIST-reproducibility
![gabst_01dic](https://github.com/user-attachments/assets/0d36443d-007f-423b-ae75-64bb1a3e23c2)
This repository includes R scripts to reproduce results, including information from figures and tables, from the original singIST paper (bioRxiv, 2024). singIST is a method for comparative single-cell transcriptomics between disease models and humans, which provides explainable and quantitative insights into single-cell transcriptomics alignment. 

## Table of Contents

1. [Repository structure](#Repository-structure)
2. [Associated code of figures and tables](#Associated-code-of-figures-and-tables)
3. [Setup Instructions](#Setup-Instructions)
4. [Requirements](#Requirements)
5. [R Dependencies](#R-Dependencies)
   

# Repository structure
The repository is organized as follows:
- `1_rawdata`: human scRNA-seq data wrangling and gene set extraction with MsigDB. Its code is in Databricks R scripts.
- `2_step1_singIST`: asmbPLS-DA model training, validity test of the optimal asmbPLS-DA, and parameter variabilities and significance.
- `3_step2_3_4_singIST`: biological link function, computation of reference recapitulation metrics, computation of predicted recapitulation metrics, and predicted recapitulation metrics as a fraction of reference recapitulations.
- `4_graphical_outputs`: graphics of superpathway recapitulation, observed one-to-one orthology, cell type recapitulation and gene contribution.
  
Each report comes with a companion folder with the exported results. 

# Associated code of figures and tables 
| Figure/Table      | R script                       |
| ------------ | --------------------------------- |
| Table 1 | asmbPLSDA_validation.R  | 
| Table 2  |  asmbPLSDA_validation.R            |
| Fig 3    |  figure_3.R             |
| Fig 4    |  figure_4.R              |
| Fig 5    |  figure_5.R               |

# Setup Instructions
1. Ask repo owner (morunoaitor@gmail.com) JSON credentials to access Google Cloud Storage raw data folder.  
2. Run following R script to download raw data:
```R
install.packages(c("Seurat", "googleCloudStorageR"))
library(googleCloudStorageR)

# Parameters
credentials_JSON = "CREDENTIAL_FILE.json" # CHANGE to JSON file provided by repo owner
file_path = "local_path/" # CHANGE to desired local path to download data
bucket_name = "human_diseasemodel_data"

# Download raw data in local folder
gcs_auth("credentials_JSON.json")
## OVA disease model
gcs_get_object("/diseasemodels/OVA.rds", bucket = bucket_name, 
               saveToDisk = paste0(file_path, "OVA.rds"), overwrite = TRUE)
## OXA+IMQ disease models
gcs_get_object("/diseasemodels/OXA_IMQ.rds", bucket = bucket_name, 
               saveToDisk = paste0(file_path, "OXA_IMQ.rds"), overwrite = TRUE)
## Human
gcs_get_object("/human/human_bangert.rds", bucket = bucket_name, 
               saveToDisk = paste0(file_path, "human_bangert.rds"), overwrite = TRUE)
```
3. Clone the repository:
```bash
git clone https://github.com/amoruno/singIST-reproducibility.git
```
# Requirements
## System Requierements
### Databricks script (.dbs)
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

### R script (.R)
```
R version 4.4.2 (2024-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 10 x64 (build 19045)

Matrix products: default


locale:
[1] LC_COLLATE=Spanish_Spain.utf8  LC_CTYPE=Spanish_Spain.utf8   
[3] LC_MONETARY=Spanish_Spain.utf8 LC_NUMERIC=C                  
[5] LC_TIME=Spanish_Spain.utf8    

time zone: Europe/Madrid
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
[1] compiler_4.4.2 tools_4.4.2   
```

# R Dependencies

Find them in DEPENDENCIES.R
