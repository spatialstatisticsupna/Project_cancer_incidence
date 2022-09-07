# Spanish Simulation Study

This folder contains the necessary R code to replicate and reproduce the Spanish Simulation Study.

In the data folder the cartography used and the Spanish data can be found. This data set consist of the variables:

 - district: The number of the region.

 - population: Population.

The Principal_SimulationStudy.Rmd file allows you to replicate the Scenarios defined in the Simulation Study.

The Original_Simulations folder contains the 100 simulations that we have used to obtain the results of the paper. With these simulations, the results of the paper can be reproduce.

The functions folder contains part of the R code necessary to implement the proposed multivariate models.

The Multivariate_Models.Rmd file allows you to fit the multivariate models proposed in the paper using INLA.

The Select_Simulations.Rmd file allows you to select 100 simulations among all simulations defined and to save only the variables of interest in our study.

The Results.Rmd file allows you to reproduce and replicate the results obtained in the paper.
