################################################################################
###############       Spatio-temporal Shared models (INLA)       ###############
################################################################################


########################################
## Data organization for INLA         ##
########################################
data <- Data_sim
results<-list()

    
## Number of areas
n<- length(unique(data$district))
    
crimes <- as.character(unique(data$disease))
J <- length(crimes)
    
##Indices:
data$ID_area<-1:(2*n)
data$alpha1 <- rep(1:0,each=n)
data$alpha2 <- rep(0:1,each=n)
data$ID_unst <- c(1:n)

    
########################################
##         d=1 lung, d=2 LOCP         ##
########################################
    
data$ID_type<-rep(c(1,2),each=n)
    
#### Definimos los PC priors
pcprec <- list(theta=list(prior='pc.prec', param=c(1,0.01)))
    
    
#### Models:
formula = observed~ -1 + alpha1 + alpha2 +
  f(ID_area, model="besag2", graph=g,hyper=pcprec, scale.model = T) + 
  f(ID_unst, model="iid",hyper=pcprec,replicate = ID_type)


lc <- c()

for(l in 1:n) {
  for (k in 1:2) {
    sh <- rep(0, 2*n)
    sh[ which(data$ID_area == l+(k-1)*n) ] <- 1
    un <- rep(0, n)
    un[ which(data$ID_unst == as.character(l))[k] ] <- 1
    lc.i <- inla.make.lincomb(ID_area = sh, ID_unst = un)
    names(lc.i) <- paste0("Lincomb.", l+(k-1)*n)
    lc <- c(lc, lc.i)
  }
}

    
result_shared_unst2 = tryCatch({
  inla(
    formula, 
    family = "poisson",
    data = data,
    E=population,
    lincomb = lc,
    control.compute=list(dic=TRUE, 
                         cpo=TRUE, 
                         waic=TRUE, 
                         hyperpar=TRUE),
    control.inla=list(strategy="simplified.laplace",
                      verbose = F, 
                      numint.maxfeval= 100000),
    control.predictor = list(link=1,compute = TRUE))
  }, error = function(msg) {
  print(msg)
  return(NULL)
  })
    




