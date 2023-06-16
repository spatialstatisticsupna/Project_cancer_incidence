# Predicting cancer incidence in regions without population-based cancer registry using mortality 

This repository contains the R code to fit in INLA the shared component models and the M-models, and the code to replicate and reproduce the simulation study described in the paper _"Predicting cancer incidence in regions without population-based cancer registry using mortality"_ [(Retegui et al., 2023)](https://doi.org/10.1093/jrsssa/qnad077).

## Table of contents

1.  [Note](#Note)
2.  [R code](#Rcode)
  1. [Germany data](#Germany)
  2. [Spanish Simulation Study](#Spain)
3.  [Acknowledgements](#Acknowledgements)
4.  [References](#Ref)

# Note <a name="Note"/>

Due to confidentiality issue the original data provided by the German Centre for Cancer Registry Data, ZfKD, used in the paper has been modified. So the code used in the paper and presented [here](https://github.com/spatialstatisticsupna/Project_cancer_incidence/tree/main/R), can be used to replicate but not to reproduce the results shown with German incidence and mortality cancer data.

# R code <a name="Rcode"/>
This folder contains the R code to replicate the spatial multivariate models and to replicate and reproduce the simulation study described in the paper.

The code to replicate the spatial multivariate models and the results obtained with German cancer data can be found [here](https://github.com/spatialstatisticsupna/Project_cancer_incidence/tree/main/R/Germany_data).

The code to replicate and reproduce the simulation study and to reproduce the results obtained in the simulation study can be found [here](https://github.com/spatialstatisticsupna/Project_cancer_incidence/tree/main/R/Spanish_Simulation_Study).

# Acknowledgements <a name="Acknowledgements"/>
We would like to thank the Centre for Cancer Registry Data (ZfKD), Germany, for providing data on federal state level.

This work has been supported by Project PID2020-113125RB-I00 (AEI), Proyecto Jóvenes Investigadores PJUPNA2018-11 and Ayudas Predoctorales Santander UPNA
2021-2022.
![plot](https://github.com/spatialstatisticsupna/Estimating_LOCP_cancer_mortality_rates/blob/main/micin-aei.jpg)

# References <a name="Ref"/>

[Retegui, G., Etxeberria, J., Riebler, A. and Ugarte, MD (2023). Predicting cancer incidence in regions without population-based cancer registry using mortality. _Journal of the Royal Statistical Society Series A (published online)_](https://doi.org/10.1093/jrsssa/qnad077)
