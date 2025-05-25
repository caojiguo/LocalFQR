This folder consists of the codes for used to perform the real data analysis in the manuscript.

###########################################################

1. soybean_fit_initial.R: the R code used to obtain initial estimate for the slope function without locally sparse feature based on the soybean yield data. It performs cross validation to decide the optimal tuning parameter value for the roughness penalty (there is no sparse penalty in this step).

2. soybean_cross_validation.R: the R code used to perform cross validation to decide the optimal values when fitting the model with both roughness penalty and group lasso type penalty. To run this file, you need to run soybean_fit_initial.R first. The running time of this file is extremely long. We usually run it in a parallel way on high performance clusters. If you want to directly reproduce the estimation results in the manuscript, you can directly run the following file. 

3. soybean_fit_sparse.R: the R code used to obtain the final estimate for the slope function with locally sparse feature based on the initial value obtained from soybean_fit_initial.R and optimal tuning parameter values decide by soybean_cross_validation.R. To run this file, you need to run soybean_fit_initial.R first.

4. nodes.RDS: contains all the nodes of the triangulation.
   
   tri.RDS: contains the information about all the triangles formulated based on the nodes information saved in "nodes.RDS".

5. soy_avg_x.RDS: contains observations on daily minimum temperature and daily maximum temperature over one year (365 days) at county level in Kansas.

   soy_avg_y.RDS: contains observations on soybean yield per unit at the end of the year at county level in Kansas.

   soy_avg_add.RDS: contains observations on annual precipitation, harvested area and irrigated area at county level in Kansas.

6. Triangulation.R: the R code used to generate nodes.RDS and tri.RDS for the triangulation. You do not need to run this since we have provided the files for triangulation.

##########################################################

To reproduce the analysis, please follow the steps below.

Steps:

A. Set working directory to this folder

B. In R, first run soybean_fit_initial.R to install necessary packages and obtain the initial estimate for the slope function.
	
C. The next step is to fit the model with both roughness and group LASSO penalties. If you want to reproduce all the cross validation results of parameter tuning for the roughness and group LASSO penalties, you should run soybean_cross_validation.R. This will take a very long time. If you want to skip the cross validation step and directly obtain the estimate of slope function using the optimal tuning parameter values, you can just run soybean_fit_sparse.R. The optimal tuning parameter values have included in soybean_fit_sparse.R.