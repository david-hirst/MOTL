# Glioblastoma use case application scripts

We downloaded the mRNA and DNA methylation zip files from [this zenodo repository](https://zenodo.org/records/7380252) in order to create a multi-omics target dataset.

The following scripts were used to carry out use case application of MOTL and comparison with direct MOFA faxtorization.

[01_TCGA_Lrn_pre_MOFA.R](https://github.com/david-hirst/MOTL/blob/main/GlioblastomaUseCase/01_TCGA_Lrn_pre_MOFA.R): This file contains the workflow for pre-processing the TCGA learning dataset before MOFA factorization

[02_Trg_pre_MOFA.R](https://github.com/david-hirst/MOTL/blob/main/GlioblastomaUseCase/02_Trg_pre_MOFA.R): This file contains the workflow for pre-processing the glioblastoma target dataset before MOFA factorization

[03_Lrn_Trg_MOFA.py](https://github.com/david-hirst/MOTL/blob/main/GlioblastomaUseCase/03_Lrn_Trg_MOFA.py): This file contains the workflow for performing MOFA factorization on the learning and target datasets

[04_TCGA_Lrn_Intercepts.R](https://github.com/david-hirst/MOTL/blob/main/GlioblastomaUseCase/04_TCGA_Lrn_Intercepts.R): This file contains the workflow for inferring intercepts for the MOFA factorization of the learning dataset

[05_Trg_post_MOFA.ipynb](https://github.com/david-hirst/MOTL/blob/main/GlioblastomaUseCase/05_Trg_post_MOFA.ipynb): This file contains the workflow for investigating the MOFA factorization of the glioblastoma target dataset

[06_Trg_MOTL.R](https://github.com/david-hirst/MOTL/blob/main/GlioblastomaUseCase/06_Trg_MOTL.R): This file contains the workflow for pre-processing the target dataset before MOTL factorization, appplication of MOTL, investigating the MOTL factorization, and comparison with the MOFA factorization

[07_orsumCommands.txt](https://github.com/david-hirst/MOTL/blob/main/GlioblastomaUseCase/07_orsumCommands.txt): This file contains the instructions for creating GSEA plots using the orsum python package
