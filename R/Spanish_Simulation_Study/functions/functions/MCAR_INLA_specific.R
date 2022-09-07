MCAR_INLA_st <- function(carto=NULL,
                         data=NULL, 
                         ID.area=NULL,
                         O=NULL,
                         E=NULL,
                         prior.spatial="Leroux",
                         prior.disease="FE",
                         model="global", 
                         strategy="simplified.laplace"
                         ){
        
        ## Check for errors ##
        if(is.null(carto))
                stop("the carto argument is missing")
        if(!any(class(carto) %in% c("SpatialPolygonsDataFrame","sf")))
                stop("the carto argument must be of class 'SpatialPolygonsDataFrame' or 'sf'")
        if(is.null(ID.area))
                stop("the ID.area argument is missing")
        
        
        if(is.null(O))
                stop("the O argument is missing")
        if(is.null(E))
                stop("the E argument is missing")
        
        if(!(prior.spatial %in% c("Leroux","intrinsic","proper","BYM","BYM2")))
                stop("invalid prior.spatial argument")
        if(!(prior.disease %in% c("FE","RE")))
                stop("invalid prior.disease argument")
        
        
        
        if(!(model %in% c("global")))  
                stop("invalid model argument")
        if(!(strategy %in% c("gaussian","simplified.laplace","laplace","adaptative")))
                stop("invalid strategy argument")
        
        
        ## Pre-processing data                  ##
        ## ---------------------------------------
        cat("STEP 1: Pre-processing data\n")
        
        
        ## Transform 'SpatialPolygonsDataFrame' object to 'sf' class ##
        carto <- sf::st_as_sf(carto)
        
        ## Order the data ##
        j <- length(data)
        if(is.null(names(data))) names(data) <- paste("disease",1:j,sep=".")
        
        for(i in names(data)){
                if(!ID.area %in% colnames(carto))
                        stop(sprintf("'%s' variable not found in carto object",ID.area))
                if(!ID.area %in% colnames(data[[i]]))
                        stop(sprintf("'%s' variable not found in '%s' data.frame object",ID.area,i))
                
                if(!O %in% colnames(data[[i]]))
                        stop(sprintf("'%s' variable not found in '%s' data.frame object",O,i))
                if(!E %in% colnames(data[[i]]))
                        stop(sprintf("'%s' variable not found in '%s' data.frame object",E,i))
                if(length(unique(data[[i]][,ID.area]))!=length(unique(sf::st_set_geometry(carto, NULL)[,ID.area])))
                        stop(sprintf("The values of '%s' variable do not match between data and carto objects", ID.area))
                if(!all(as.character(sort(unique(data[[i]][,ID.area])))==as.character(sort(unique(sf::st_set_geometry(carto, NULL)[,ID.area])))))
                        stop(sprintf("The values of '%s' variable do not match between data and carto objects", ID.area))
        }
        
        ## Define hyperprior distributions      ##
        ## ---------------------------------------
        sdunif="expression:
          logdens=-log_precision/2;
          return(logdens)"
        
        lunif = "expression:
          a = 1;
          b = 1;
          beta = exp(theta)/(1+exp(theta));
          logdens = lgamma(a+b)-lgamma(a)-lgamma(b)+(a-1)*log(beta)+(b-1)*log(1-beta);
          log_jacobian = log(beta*(1-beta));
          return(logdens+log_jacobian)"
        
        
        ## Formula for INLA model               ##
        ## ---------------------------------------
        form <- "O ~ -1+"
        form <- paste(form, paste(paste0("I",1:j),collapse="+"),sep="")
        form <- paste(form, "+ f(idx, model=Mmodel, constr=FALSE, extraconstr=list(A=A.constr.s, e=rep(0,j)))")
        if(prior.spatial %in% c("BYM")){
                form <- paste(form, "+ f(idx.v, model=Mmodel.v, constr=FALSE, extraconstr=list(A=A.constr.s, e=rep(0,j)))")
        }
           
        formula <- stats::as.formula(form)
        
        
        
        ## Global model                         ##
        ## ---------------------------------------
        if(model=="global"){
                
                ## Fitting global model with INLA       ##
                ## ---------------------------------------
                cat("STEP 2: Fitting model with INLA (this may take a while...)\n")
                
                ## order
                data <- lapply(data, function(x) x[order(x[,ID.area]),] )
                
                n <- length(unique(data[[1]][,ID.area]))
                
                
                ## Adjacency matrices                   ##
                ## ---------------------------------------
                
                ## W: adjacency matrix (spatial)
                adj <- spdep::poly2nb(carto)
                Ws <- as(nb2mat(adj, style = "B"), "Matrix")
                
                
                ## Precision matrices                   ##
                ## ---------------------------------------
                Rs <- Matrix::Diagonal(n,colSums(Ws))-Ws
                Rs.Leroux <- Matrix::Diagonal(n)-Rs
                
                ## Define appropriate constraints matrices##
                ## ---------------------------------------
                A.constr.s<- kronecker(diag(j), matrix(1,1,n))
                
                
                ## Data for INLA model                  ##
                ## ---------------------------------------
                data.INLA <- data.frame(O=as.numeric(unlist(lapply(data, function(x) x[,O]))),
                                        E=as.numeric(unlist(lapply(data, function(x) x[,E]))),
                                        Area=as.character(unlist(lapply(data, function(x) x[,ID.area]))),
                                        ID.area=1:n,
                                        ID.disease=rep(1:j, each=n)
                                        )
                intercepts <- fastDummies::dummy_cols(data.INLA$ID.disease)[,-1]
                intercepts[intercepts==0] <- NA
                colnames(intercepts) <- paste0("I",1:j)
                data.INLA <- cbind(data.INLA, intercepts)
                
                data.INLA$idx <- (data.INLA$ID.disease-1)*n + data.INLA$ID.area
                data.INLA$idx.v <- data.INLA$idx
                
                
                ## Initial values for hyperparameters   ##
                ## ---------------------------------------
                library(tidyverse)
                
                ## Spatial initial values (hyperparameters)
                sir.s <- data.INLA %>% 
                        dplyr::group_by(ID.disease, ID.area) %>%
                        dplyr::summarise(SIR = sum(O)/sum(E))
                Sigma.s <- cov(matrix(sir.s$SIR, nrow=n, ncol=j, byrow = FALSE), use="complete.obs")
                Rho.s <- cor(matrix(sir.s$SIR, nrow=n, ncol=j, byrow = FALSE), use="complete.obs")
                
                
                initial.values.s <- as.vector(
                        c(-log(diag(Sigma.s)),
                          log(1+Rho.s[lower.tri(Rho.s, diag = FALSE)])-log(1-Rho.s[lower.tri(Rho.s, diag = FALSE)])
                          ))
                
                
                ## Define selected Mmodel               ##
                ## ---------------------------------------
                ## Spatial
                if(prior.spatial=="intrinsic" & prior.disease=="FE"){
                        Mmodel <- INLA::inla.rgeneric.define(Mmodel_icar_fe, debug=FALSE, k=j, W=Ws, initial.values=initial.values.s)
                }
                if(prior.spatial=="intrinsic" & prior.disease=="RE"){
                        Mmodel <- INLA::inla.rgeneric.define(Mmodel_icar_re, debug=FALSE, k=j, W=Ws, initial.values=initial.values.s)
                }
                if(prior.spatial=="Leroux" & prior.disease=="FE"){
                        Mmodel <- INLA::inla.rgeneric.define(Mmodel_lcar_fe, debug=FALSE, k=j, W=Ws, initial.values=initial.values.s, alpha.min=0, alpha.max=1)
                }
                if(prior.spatial=="Leroux" & prior.disease=="RE"){
                        Mmodel <- INLA::inla.rgeneric.define(Mmodel_lcar_re, debug=FALSE, k=j, W=Ws, initial.values=initial.values.s, alpha.min=0, alpha.max=1)
                }
                if(prior.spatial=="proper" & prior.disease=="FE"){
                        Mmodel <- INLA::inla.rgeneric.define(Mmodel_pcar_fe, debug=FALSE, k=j, W=Ws, initial.values=initial.values.s, alpha.min=0, alpha.max=1)
                }
                if(prior.spatial=="proper" & prior.disease=="RE"){
                        Mmodel <- INLA::inla.rgeneric.define(Mmodel_pcar_re, debug=FALSE, k=j, W=Ws, initial.values=initial.values.s, alpha.min=0, alpha.max=1)
                }
                if(prior.spatial=="BYM" & prior.disease=="FE"){
                        Mmodel <- INLA::inla.rgeneric.define(Mmodel_icar_fe,  debug=FALSE, k=j, W=Ws, initial.values=initial.values.s)
                        Mmodel.v <- INLA::inla.rgeneric.define(Mmodel_iid_fe, debug=FALSE, k=j, W=Ws, initial.values=initial.values.s)
                }
                if(prior.spatial=="BYM" & prior.disease=="RE"){
                        Mmodel <- INLA::inla.rgeneric.define(Mmodel_icar_re, debug=FALSE, k=j, W=Ws, initial.values=initial.values.s)
                        Mmodel.v <- INLA::inla.rgeneric.define(Mmodel_iid_re, debug=FALSE, k=j, W=Ws, initial.values=initial.values.s)
                }
                if(prior.spatial=="BYM2" & prior.disease=="FE"){
                        Mmodel <- INLA::inla.rgeneric.define(Mmodel_bym2_fe, debug=FALSE, k=j, W=W, initial.values=initial.values.s, alpha.min=0, alpha.max=1)
                }
                if(prior.spatial=="BYM2" & prior.disease=="RE"){
                        Mmodel <- INLA::inla.rgeneric.define(Mmodel_bym2_re, debug=FALSE, k=j, W=W, initial.values=initial.values.s, alpha.min=0, alpha.max=1)
                }
                
                
                ## Fit the INLA model                   ##
                ## ---------------------------------------
                Model <- INLA::inla(formula, family="poisson", data=data.INLA, E=E,
                                    control.predictor=list(link=1, compute=TRUE, cdf=c(log(1))),
                                    control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, config=TRUE),
                                    control.inla=list(strategy=strategy),
                                    debug = F,
                                    verbose = F)
                
                Model$Mmodel <- list(model=model, prior.spatial=prior.spatial, prior.disease=prior.disease)
                
        }
        return(Model)
}
