
# Simulation study scripts

These are the scripts for reproducing the simulation study, listed in the order in which they are run.

[Multi_Omics_Simulator.py](https://github.com/david-hirst/MOTL/blob/main/SimulationStudy/Multi_Omics_Simulator.py): This contains the function for generating a multi-omics dataset based on a groundtruth factor structure

[01_Simulation_Loop.py](https://github.com/david-hirst/MOTL/blob/main/SimulationStudy/01_Simulation_Loop.py): This imports Multi_Omics_Simulator.py, and then simulates a specified number of multi-omics datasets for a given simulation configuration. The configurations and the number of simulated datasets is set by the user

[02_FeatureVarsEval.Rmd](https://github.com/david-hirst/MOTL/blob/main/SimulationStudy/02_FeatureVarsEval.Rmd):This loops through each simulated multi-omics data set for a given simulation configuration. For each dataset it checks the featurewise variance and outputs a file with the variance for each feature, as well as a file with rankings based on the variance. It does this for the learning dataset and the target dataset.

[03_FactorVarsEval.Rmd](https://github.com/david-hirst/MOTL/blob/main/SimulationStudy/03_FactorVarsEval.Rmd): This loops though each simulation instance for a given configuration. It calculates the number of active factors underlying learning and target datasets and outputs a fies with the active factor count by simulation instance number

[04_FctrznMOFA_ActiveK.py](https://github.com/david-hirst/MOTL/blob/main/SimulationStudy/04_FctrznMOFA_ActiveK.py): This loops through each simulation instance of a given configuration and factorizes either the learning or target dataset with MOFA. The user specifies which dataset to factorize. The factorization result is ouputted to a named folder

[05_Lrn_Intcpt_MLE.Rmd](https://github.com/david-hirst/MOTL/blob/main/SimulationStudy/05_Lrn_Intcpt_MLE.Rmd): This loops through simulation configurations and simulation instances. For each learning dataset it estimates an intercept for the corresponding MOFA factorization so that an intercept can be used in transfer learning.

[06_TL_VI.Rmd](https://github.com/david-hirst/MOTL/blob/main/SimulationStudy/06_TL_VI.Rmd): This loops through simulation configurations and instances and performs transfer learning with variational inference

[06b_TL_VI_permute.Rmd](https://github.com/david-hirst/MOTL/blob/main/SimulationStudy/06b_TL_VI_permute.Rmd): This loops through simulation configurations and instances and performs transfer learning with variational inference after permuting the values in proportions of the vectors in the W matrix inferred from L

[06c_intNMF.Rmd](https://github.com/david-hirst/MOTL/blob/main/SimulationStudy/06c_intNMF.Rmd): This loops through simulation configurations and instances and factorizes the target dataset with IntNMF

[06d_MoCluster.Rmd](https://github.com/david-hirst/MOTL/blob/main/SimulationStudy/06d_MoCluster.Rmd): This loops through simulation configurations and instances and factorizes the target dataset with MoCluster

[07_FctrznEvaluation.Rmd](https://github.com/david-hirst/MOTL/blob/main/SimulationStudy/07_FctrznEvaluation.Rmd): This loops through each simulation configuration, instance and factorization and calculates some evaluation metrics. An rds file is outputted for each simulation configuration

[08_FctrznEvaluationCmbnd.Rmd](https://github.com/david-hirst/MOTL/blob/main/SimulationStudy/08_FctrznEvaluationCmbnd.Rmd): This loops through simulation configurations and combines all the factorization evaluation scores into a single data frame. This is used for plotting and chekcing statistical significance
