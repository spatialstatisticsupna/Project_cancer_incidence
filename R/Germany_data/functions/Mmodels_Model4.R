################################################################################

library(tidyverse)

## Packages loading
packages <- c("INLA", "stats", "spdep","sf", "fastDummies")
invisible(lapply(packages, library, character.only = TRUE))


################################################################################
## Data organization for INLA                                                 ##
################################################################################
## Load carto
Carto_UP <-  carto_use

################################################################################
##  1) Fitting models                                                         ##
################################################################################
source("./functions/functions/Mmodel_icar.R") ## Intrinsic multivariate CAR latent effect
source("./functions/functions/Mmodel_lcar.R") ## Leroux multivariate CAR latent effect
source("./functions/functions/Mmodel_pcar.R") ## Proper multivariate CAR latent effect
source("./functions/functions/Mmodel_iid.R")  ## iid
source("./functions/functions/Mmodel_bym2.R") ## BYM2

## Fit a spatial-temporal multivariate Poisson mixed model to areal count data, 
## where dependence between spatial/temporal patterns of the diseases is addressed
## through the use of M-models
source("./functions/functions/MCAR_INLA_specific.R")


## The strategy to use for the approximations
strategy <- "simplified.laplace"

priors <- c("iCAR","LCAR","pCAR","BYM","BYM2")
method <- c("FE","RE")

prior_spatial <- c("intrinsic","Leroux","proper","BYM","BYM2")

Data_UP <- Data_sim2

for (p in 1:length(prior_spatial)) {
  for (m in 1:length(method)) {
    eval(parse(text = paste0(
      "result_",priors[p],"_",method[m],"_specI <- tryCatch({
  MCAR_INLA_st(
    carto = Carto_UP,
    data = Data_UP,
    ID.area = 'ID_1',
    O = 'counts',
    E = 'population',
    prior.spatial = prior_spatial[p],
    prior.disease = method[m],
    model = 'global',
    strategy = strategy
  )
}, error = function(msg) {
      print(msg)
      return(NULL)
})"
    )))
  }
}
    
