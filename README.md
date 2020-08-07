# Neurodegenerative Disease Progression Model for Parkinson's Disease Analysis
Yuan Zhou

## Introduction

Parkinson's disease (PD) progression can be studied using DaTscans collected from living PD subjects, which reflect the density of dopamine transporters in the striatum in the brain, hence measure the availability of dopaminergic neurons left in the PD subjects.

The PPMI dataset provides 449 subjects that have DaTscans collected at baseline, year 1, 2, 4, and 5. By applying a mixture of linear dynamical system (MLDS) to this longitudinal dataset, we can analyze the progression patterns and subtypes of PD.

This repository contains the Matlab implementation of this MLDS algorithm. 

If you find some of the code helpful, please cite the corresponding papers.  

Zhou, Y., Tagare, H. D., "Bayesian Longitudinal Modeling of Early Stage Parkinson’s Disease Using DaTscan Images", *International Conference on Information Processing in Medical Imaging*, Hongkong, 2020  
(https://link.springer.com/chapter/10.1007/978-3-030-20351-1_31)  

Zhou, Y., Tinaz, S., Tagare, H. D., "Robust Bayesian Analysis of Early-Stage Parkinson's Disease Progression Using DaTscan Images". *IEEE Transactions on Medical Imaging*, 2020, Under Review  

#### Organization

The folder "mlds" contains the algorithm implementation of MLDS.  

The folder "utils" contains some utility functions.  

## Mixture of Linear Dynamical System

The MLDS models the exponential decay of signals in the DaTscans using a linear dynamical system (LDS) with a centrosymmetric constraint on the transition matrix. Extending this LDS to K different transition matrices, we can capture the subtypes inherent in the dataset.
  

![mlds-flowchart](./figure/overview_mlds.jpg)

To see the demo, change the directory to "mlds" in Matlab and run
```
run_algo;
```


#### Organization

The "EM" folder contains the EM algorithm to solve the MLDS.

The "gibbs_sampling" folder contains the Gibbs sampling algorithm to solve the MLDS.

The "model_selection" folder contains the model selection (Bayesian or cross validation) code for both the Gaussian noise case and the *t*-distribution noise case.

The "utility" folder contains the utility code used for the MLDS.


## Contact

If you have any questions, please contact:

Yuan Zhou  
Department of Radiology and Biomedical Imaging  
Yale School of Medicine  

zhouyuanzxcv@gmail.com

