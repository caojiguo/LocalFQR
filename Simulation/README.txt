This folder consists of the codes for the algorithm used to perform the simulation studies in the manuscript.

You will need the following R files to run the simulations. In this README file, I will explain each R file and then outline how to run them in the correct order.

###########################################################

1. DataGen_I.R - DataGen_VII.R: the R code used to generate data for Scenario I-VII.
   
2. sfqr_fit_initial.R: the R code used to obtain the initial estimates of the parameters, which fits the model without sparsity penalties.

3. sfqr_fit_sparse.R: the R code used to first identify the signal regions of $\beta(t,u)$ using the proposed locally sparse estimation method with both roughness and group lasso type penalties based any given tuning parameter values and the initial estimates obtained from sfqr_fit_initial.R, and then refit the model again only using the data within the identified signal regions. Usually we need to conduct cross validation to decide the optimal tuning parameter values, but in this file, the tuning parameters are manually input by the user. To run this file, you need to run sfqr_fit_initial.R first.

4. nodes.RDS: contains all the nodes of the triangulation.
   
   tri.RDS: contains the information about all the triangles formulated based on the nodes information saved in "nodes.RDS".


5. Triangulation.R: the R code used to generate nodes.RDS and tri.RDS for the triangulation. You do not need to run this since we have provided the files for triangulation.

##########################################################

To reproduce the analysis, please follow the steps below.

Steps:

A. Set the working directory to this folder.

B. In R, first run sfqr_fit_initial.R to install necessary packages, generate simulated data and obtain the initial estimate for the slope function.
	
C. The next step is to fit the model with both roughness and group LASSO penalties. Run sfqr_fit_sparse.R to obtain the locally sparse estimate for the slope function for any given tuning parameter values.