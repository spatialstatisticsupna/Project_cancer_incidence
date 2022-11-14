# German data

This folder contains the necessary R code to fit the multivariate models used in the paper and the results obtained. Due to confidentiality issues, incidence and mortality data from Germany are not shown in the repository. Therefore, some modifications have been made to the data in order to replicate the data analysis.

In the data folder one may find the cartography and the modified German incidence and mortality data for each cancer location: lung and breast cancer (females), and prostate cancer. Files named  lungF.rda,  breast.rda, and prostate.rda, respectively. The data files have the following varibales:

 - ID_area: a region index. Takes values 1 to 16

 - population: population

 - counts: observed incidence and mortality counts (one after the other)

- health_outcome: takes value "inc" for cancer incidence and value "mort" for cancer mortality.

Part of the R code needed to implement the proposed multivariate models can be found in the functions folder.

The Multivariate_Models.Rmd file allows you to fit the multivariate models proposed in the paper using INLA.

The Results.Rmd file allows you to reproduce the results obtained in the paper.
