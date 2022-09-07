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
data$ID_unst1 <- c(rep("NA",n),1:n)
data$ID_unst2 <- c(1:n, rep("NA",n))
    
########################################
##         d=1 lung, d=2 LOCP         ##
########################################
    
data$ID_type<-rep(c(1,2),each=n)
    
#### Definimos los PC priors
pcprec <- list(theta=list(prior='pc.prec', param=c(1,0.01)))
    
    
#### Models:
formula = observed~ -1 + alpha1 + alpha2 +
  f(ID_area, model="besag2", graph=g,hyper=pcprec, scale.model = T) + 
  f(ID_unst1, model="iid",hyper=pcprec) +
  f(ID_unst2, model="iid",hyper=pcprec)

# lc <- c()
# 
# for(l in 1:n) {
#   sh <- rep(0, 2*n) #vector de 0s para cada categoria
#   sh[ which(data$ID_area == l+n) ] <- 1 #que categoria queremos tener en cuenta y con que factor (1)
#   un <- rep(0, n+1)
#   factor.un <- levels(factor(data$ID_unst))
#   un[ which(factor.un == as.character(l)) ] <- 1
#   lc.i <- inla.make.lincomb(ID_area = sh, ID_unst = un)
#   names(lc.i) <- paste0("Lincomb.", l)
#   lc <- c(lc, lc.i)
# }

    
result_shared_unst2_indep = tryCatch({
  inla(
    formula, 
    family = "poisson",
    data = data,
    E=population,
    # lincomb = lc,
    control.compute=list(dic=TRUE, 
                         cpo=TRUE, 
                         waic=TRUE, 
                         hyperpar=TRUE,
                         config = TRUE),
    control.inla=list(strategy="simplified.laplace",
                      verbose = F, 
                      numint.maxfeval= 100000),
    control.predictor = list(link=1,compute = TRUE))
  }, error = function(msg) {
  print(msg)
  return(NULL)
  })
    




