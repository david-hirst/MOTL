# TCGA evaluation study scripts

This folder contains the functions and workflows used to carry out the TCGA evaluation study.

## Workflows for data acquisition, pre-processing and factorization
[TCGA_downloadData.R](https://github.com/david-hirst/MOTL/blob/main/TCGAStudy/TCGA_downloadData.R): This file contains the workflow for downloading and saving TCGA multi-omics data

[TCGA_preprocessedData.R](https://github.com/david-hirst/MOTL/blob/main/TCGAStudy/TCGA_preprocessedData.R): This file contains examples of workflows for creating and pre-processing TCGA Learning, reference, and target datasets

[TCGA_MOFA.py](https://github.com/david-hirst/MOTL/blob/main/TCGAStudy/TCGA_MOFA.py): This file contains workflows for factorizating TCGA learning, reference and target datasets with MOFA

[TCGA_Lrn_Intercepts.R](https://github.com/david-hirst/MOTL/blob/main/TCGAStudy/TCGA_Lrn_Intercepts.R): This file contains the workflow for inferring intercepts for a MOFA factorization of a TCGA learning dataset

[TL_VI.R](https://github.com/david-hirst/MOTL/blob/main/TCGAStudy/TL_VI.R): This file contains the worflow for applying MOTL to a TCGA target dataset to carry out matrix factorization with transfer learning

## Workflows for evaluation measures and figures
NOTE: THE CODE IN THE PREVIOUS SECTION IS AN OPTIMIZED VERSION OF THE CODE WE USED TO CARRY OUT THE TCGA EVALUATION STUDY. THE CODE WE USED FOR EVALUATING THE FACTORIZATIONS IS STILL UNDERGOING THIS PROCESS.

[TCGA_FctrznEvaluation.R](https://github.com/david-hirst/MOTL/blob/main/TCGAStudy/TCGA_FctrznEvaluation.R): This file contains the workflow used for calulating evaluation measures for factorizations of TCGA target datasets.

[TCGA_FctrznEvaluation_Summary.R](https://github.com/david-hirst/MOTL/blob/main/TCGAStudy/TCGA_FctrznEvaluation_Summary.R): This file contains the workflow for plotting and performing statistical testing on the TCGA evaluation measure values

[TCGA_FactorMatching.ipynb](https://github.com/david-hirst/MOTL/blob/main/TCGAStudy/TCGA_FactorMatching.ipynb): This file contains the workflow we used for associating differentially active groundtruth TCGA factors with genesets.

