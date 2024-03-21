# Glioblastoma use case application scripts

We downloaded the mRNA and DNA methylation zip files from [this zenodo repository](https://zenodo.org/records/7380252) in order to create a multi-omics target dataset.

The following scripts were used to carry out use case application of MOTL and comparison with direct MOFA faxtorization.

[TCGA_Lrn_pre_MOFA.R](https://github.com/david-hirst/MOTL/blob/main/GlioblastomaUseCase/TCGA_Lrn_pre_MOFA.R): This file contains the workflow for pre-processing the TCGA learning dataset before MOFA factorization

[Trg_pre_MOFA.R](https://github.com/david-hirst/MOTL/blob/main/GlioblastomaUseCase/Trg_pre_MOFA.R): This file contains the workflow for pre-processing the glioblastoma target dataset before MOFA factorization

[Lrn_Trg_MOFA.py](https://github.com/david-hirst/MOTL/blob/main/GlioblastomaUseCase/Lrn_Trg_MOFA.py): THis file contains functions ad workflow for performing MOFA factorization on the learning and target datasets

[TCGA_Lrn_Intercepts.R](https://github.com/david-hirst/MOTL/blob/main/GlioblastomaUseCase/TCGA_Lrn_Intercepts.R): This file contains the workflow for inferring intercepts for the MOFA factorization of the learning dataset

[Trg_post_MOFA.ipynb](https://github.com/david-hirst/MOTL/blob/main/GlioblastomaUseCase/Trg_post_MOFA.ipynb): This file contains the workflow for investigating the MOFA factorization of the glioblastoma target dataset

[Trg_MOTL.R](https://github.com/david-hirst/MOTL/blob/main/GlioblastomaUseCase/Trg_MOTL.R): This file contains the workflow for pre-processing the target dataset before MOTL factorization, appplication of MOTL, investigating the MOTL factorization, and comparison with the MOFA factorization

