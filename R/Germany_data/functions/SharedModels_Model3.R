################################################################################
###############       Spatio-temporal Shared models (INLA)       ###############
################################################################################


########################################
## Data organization for INLA         ##
########################################
results<-list()

Data_sim$ID_area <- (1:(2*n))
Data_sim$alpha1 <- rep(1:0,each=n)
Data_sim$alpha2 <- rep(0:1,each=n)
Data_sim$ID_unst1 <- c(1:n,rep(NA,n))
Data_sim$ID_unst2 <- c(rep(NA,n),1:n)

    
########################################
##         d=1 lung, d=2 LOCP         ##
########################################
#### Definimos los PC priors
pcprec <- list(theta=list(prior='pc.prec', param=c(1,0.01)))
    
    
#### Models:
formula = counts~ -1 + alpha1 + alpha2 +
  f(ID_area, model="besag2", graph=g,hyper=pcprec, scale.model = T) + 
  f(ID_unst1, model="iid",hyper=pcprec) + 
  f(ID_unst2, model="iid",hyper=pcprec)


# lc <- c()

# for(l in 1:n) {
#   for (k in 1:2) {
#     sh <- rep(0, 2*n)
#     sh[ which(Data_sim$ID_area == l+(k-1)*n) ] <- 1
#     un <- rep(0, n)
#     un[ which(Data_sim$ID_unst == as.character(l))[k] ] <- 1
#     lc.i <- inla.make.lincomb(ID_area = sh, ID_unst = un)
#     names(lc.i) <- paste0("Lincomb.", l+(k-1)*n)
#     lc <- c(lc, lc.i)
#   }
# }

    
result_shared_unst2_indep = tryCatch({
  inla(
    formula, 
    family = "poisson",
    data = Data_sim,
    E=population,
    # lincomb = lc,
    control.compute=list(dic=TRUE, 
                         cpo=TRUE, 
                         waic=TRUE, 
                         hyperpar=TRUE,
                         config = TRUE),
    control.inla=list(strategy="simplified.laplace",
                      int.strategy="grid",
                      verbose = F, 
                      numint.maxfeval= 100000),
    control.predictor = list(link=1,compute = TRUE))
  }, error = function(msg) {
  print(msg)
  return(NULL)
  })
    




