This folder consists of the codes for the algorithm used to perform the simulation studies in the manuscript.

You will need the following R files to run the simulations. In this README file, I will explain each R file and then outline how to run them in the correct order.

###########################################################

0. nodes.RDS: contains all the nodes of the triangulation used in the simulations.
   tri.RDS: contains the information about all the triangles formulated based on the nodes information saved in "nodes.RDS".

1. DataGen_simple_1.R: the R code used to generate data for Scenario 1 of the simulation study presented in the manuscript.
   
   DataGen_simple_2.R: the R code used to generate data for Scenario 2 of the simulation study presented in the manuscript.

   DataGen_complex.R: the R code used to generate data for Scenario 3 of the simulation study presented in the manuscript.

2. fit_sfqr_step1.R: the R code used to obtain the initial estimates of the parameters, which fits the model without sparsity penalties. (The cross validation procedure for tuning parameter is also included in this file)

3. fit_sfqr_step2.R: the R code used to first identify the inactive regions of $\beta(t,u)$ using group lasso type penalties based the initial estimates of the parameters obtained from fit_sfqr_step1.R, and then the model is refitted only using the information within the active regions. (The cross validation procedure for tuning parameters is also included in this file)


##########################################################

To reproduce the analysis, please follow the steps below.

Note that the following steps only present how to obtain one replicate of the simulation (simulate one time).

Regarding the simulation results in the manuscript, we repeat the following steps for 100 times under each simulation setting by using high performance computing clusters.

We recommend to perform the following steps on high performance computing clusters due to two reasons:

  A. The computational cost for each replication on a standard 2019 MacBook Air is more than one hour because we need to tune multiple tuning parameters using cross validation.

  B. On the clusters, we can use parallel computing to perform the simulations.

Steps:

A. In R, first install all the necessary packages. Please note that we use the package CVXR and solver MOSEK to solve optimization problem that involves group-lasso type penalties. The installation guide for MOSEK can be found at https://docs.mosek.com/latest/install/installation.html


B. Open the R files "fit_sfqr_step1.R" and "fit_sfqr_step2.R", and then modify the data generating model (DataGen_simple_1.R, DataGen_simple_2.R or DataGen_complex.R) and sample size.

C. In R, set the current folder as the working directory, and issue

   ```
   source('fit_sfqr_step1.R')

   ```

D. Next issue
 
   ```
   source('fit_sfqr_step2.R')

   ```