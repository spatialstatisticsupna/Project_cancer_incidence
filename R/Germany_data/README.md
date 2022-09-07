# German data

This folder contains the necessary R code to fit the multivariate models used in the paper and the results obtained. Due to confidentiality, incidence and mortality data from Germany are not shown in the repository. Therefore, some modifications have been made to the data in order to replicate the analysis performed.

In the Data folder the cartography used and the modified German incidence and mortality data for each cancer location analysed can be found. The cancer locations analysed have been lung cancer in females, prostate cancer and breast cancer, for which data can be found in lungF.rda, prostate.rda and breast.rda respectively. These data sets consist of the variables:

 - ID_area: a region index. Takes values 1 to 16.

 - population: population

- counts: observed incidence and mortality counts.

- health_outcome: takes value "inc" for cancer incidence and value "mort" for cancer mortality.

Part of the R code needed to implement the proposed multivariate models can be found in the functions folder.

The Multivariate_Models.Rmd file allows you to fit the multivariate models proposed in the paper using INLA.

The Results.Rmd file allows you to reproduce the results obtained in the paper.
