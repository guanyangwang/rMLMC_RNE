### Overview

The .R files contained in this folder contain the code for reproducing the numerical results in the paper 'Optimal randomized multilevel Monte Carlo for repeatedly nested expectations'

Our numerical results are provided in Section 3 (Figure 1-2), Appendix E (Table 2) and Appendix F (Figure 3-5). These experiments compare two variants of NMC estimators and READ estimator proposed in our paper in two examples. 

### File Explanation 

+ Files 'mlmc.R' and 'vanilla_NMC.r' contain the core functions of the recursive MLMC method (proposed in our paper) and the NMC method (analyzed in Rainforth et al.) respectively.

+ Files 'sin_gaussian_setup.R', 'sin_t_setup.R', and 'sigmoid_gaussian_setup.R' contain the setup of the three experiments specified in Section 3 and Appendix F of our paper.

+ Files 'sin_gaussian_exp.R', 'sin_t_exp.R', and 'sigmoid_gaussian_exp.R' contain executable R scripts to reproduce the simulation results.

### Usage

One can run 'sin_gaussian_exp.R' to reproduce the results in Section 3.1 and Appendix E, including Figure 1 and Table 1. One can run 'sin_t_exp.R' to reproduce the results in Section 3.2. One can run 'sigmoid_gaussian_exp.R' to reproduce the results in Appendix F, including Figure 3-5.  
