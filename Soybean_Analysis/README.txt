This folder consists of the codes for used to perform the real data analysis in the manuscript.

###########################################################

0. nodes.RDS: contains all the nodes of the triangulation used in the simulations.
   tri.RDS: contains the information about all the triangles formulated based on the nodes information saved in "nodes.RDS".

1. soy_x.RDS: contains observations on daily minimum temperature and daily maximum temperature over one year (365 days) at county level in Kansas.
   soy_y.RDS: contains observations on soybean yield per unit at the end of the year at county level in Kansas.
   soy_add.RDS: contains observations on annual precipitation, harvested area and irrigated area at county level in Kansas.

2. cv_output.RDS: contains the cross validation results for tuning parameter selection.

3. fit_soybean.R: the R code used to apply the proposed method to the soybean yield data, and generate the plots presented in the manuscript


##########################################################

To reproduce the analysis, please follow the steps below.

Steps:

A. In R, first install all the necessary packages. Please note that we use the package CVXR and solver MOSEK to solve optimization problem that involves group-lasso type penalties. The installation guide for MOSEK can be found at https://docs.mosek.com/latest/install/installation.html

B. In R, set the current folder as the working directory, and issue

   ```
   source('soybean_fit.R')

   ```