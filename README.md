# singIST-reproducibility
![gabst_01dic](https://github.com/user-attachments/assets/0d36443d-007f-423b-ae75-64bb1a3e23c2)
This repository includes R scripts to reproduce all steps of singIST, including information from figures and tables from the original paper (bioRxiv, 2024). singIST is a method for comparative single-cell transcriptomics between disease models and humans, which provides explainable and quantiative insights into single-cell transcriptomics alignment. 

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

