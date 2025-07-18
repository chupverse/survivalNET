predict.survivalNET <- function(object, type="net", newdata=NULL, newtimes=NULL,
                                ratetable = NULL, method = NULL, ...){
  
  if(!(type %in% c("net","lp","overall")))  stop("Argument 
                  'type' must be 'net', 'lp' or 'overall' ")
  if(type == "overall" && "m" %in% names(object))stop("The 'overall' survival prediction
  for survivalFLEXNET is still under development. Please
                            use 'type = 'net'. ")
  if(type == "lp" && "m" %in% names(object))stop("The 'lp' survival prediction
  for survivalFLEXNET is still under development. Please
                            use 'type = 'net'. ")
  ###pour retrouver progress overall -> bacupGITHUB 28_04
  if(is.null(newtimes))  { 
    newtimes <- 1:max(object$y[,1])
    
    newtimesSave <- NULL
    
  }else{
    
    newtimesSave <- newtimes
    
    newtimes <- unique( sort( c( (1:max(object$y[,1])),unique(newtimes) ) ) )
    }
  
  if(0 %in% newtimes){
    newtimes <- sort(newtimes[-(newtimes == 0)])
  }
  
  if(!is.null(newdata))
  { 
    if(!is.data.frame(newdata)) stop("Argument 'newdata' must be a data frame")
    n <- dim(newdata)[1]
    
    covnames <- names(as.data.frame(object$x))
    if ("xlevels" %in% names(object)) {
      covnames <- c(covnames, names(object$xlevels)) 
    }
    
    indic <- covnames %in% names(newdata) 
    if( sum(!indic) > 0 ) stop("Missing predictor in the data frame")
    covariates <- newdata[,covnames, drop =FALSE]
    names(covariates) <- covnames
    
  }
  
  if(is.null(newdata))  { 
    
    n <- dim(object$y)[1]
    
    covnames <- names(as.data.frame(object$x))
    covariates <- object$x 
    if ("xlevels" %in% names(object)) {
      covnames <- c(covnames, names(object$xlevels)) 
      covariates <- cbind(covariates, object$levelsval)
      colnames(covariates)[ncol(covariates)] <- names(object$xlevels)
      covariates <- as.data.frame(covariates)
    }
    
  }
  
  ### regression coefficients
  
  if ("dist" %in% names(object)) {
    
    if (!("xlevels" %in% names(object))) {   
      
      if(object$dist=="genweibull")  {
        beta <- unname( object$coefficients[1:(dim(object$x)[2])] )
        sigma <- unname( exp(object$coefficients[(dim(object$x)[2])+1]) )
        nu <- unname( exp(object$coefficients[(dim(object$x)[2])+2]) )
        theta <- unname( exp(object$coefficients[(dim(object$x)[2])+3]) )
      }
      
      if(object$dist=="weibull")  {
        beta <- unname( object$coefficients[1:(dim(object$x)[2])] )
        sigma <- unname( exp(object$coefficients[(dim(object$x)[2])+1]) )
        nu <- unname( exp(object$coefficients[(dim(object$x)[2])+2]) )
        theta <- 1
      }
      
      if(object$dist=="exponential")  {
        beta <- unname( object$coefficients[1:(dim(object$x)[2])] )
        sigma <- unname( exp(object$coefficients[(dim(object$x)[2])+1]) )
        nu <- 1
        theta <- 1
      }
    }
    if ("xlevels" %in% names(object)) {   
      
      if(object$dist=="genweibull")  {
        beta <- unname( object$coefficients[1:(dim(object$x)[2])] )
        sigmas <- unname( exp(object$coefficients[grep("^log sigma_", 
                                                       names(object$coefficients))]) ) 
        nus <- unname( exp(object$coefficients[grep("^log nu_", 
                                                    names(object$coefficients))]) )
        thetas <- unname( exp(object$coefficients[grep("^log theta_", 
                                                       names(object$coefficients))]) )
      }
      if(object$dist=="weibull")  {
        beta <- unname( object$coefficients[1:(dim(object$x)[2])] )
        sigmas <- unname( exp(object$coefficients[grep("^log sigma_", 
                                                       names(object$coefficients))]) )
        nus <- unname( exp(object$coefficients[grep("^log nu_", 
                                                    names(object$coefficients))]) )
        thetas <- rep(1,length(sigmas))
      }
      if(object$dist=="exponential")  {
        beta <- unname( object$coefficients[1:(dim(object$x)[2])] )
        sigmas <- unname( exp(object$coefficients[grep("^log sigma_", 
                                                       names(object$coefficients))]) )
        nus <- rep(1,length(sigmas))
        thetas <- rep(1,length(sigmas))
      }
    }
  }
  
  if("m" %in% names(object)){
    
    if(!("xlevels" %in% names(object))){
      m <- object$m
      mpos <- object$mpos
      beta <- unname( object$coefficients[1:(dim(object$x)[2])] )
      gamma <- unname( object$coefficients[(dim(object$x)[2]+1):(dim(object$t.table)[1])] ) 
     
    }
    if("xlevels" %in% names(object)){
      m <- object$m
      mpos <- object$mpos
      beta <- c()
      gamma <- c()
      for (i in 1:length(object$correstab)) {
        beta <- c(beta,unname( object$coefficients[(1:(dim(object$x)[2])+(dim(object$x)[2]+m+2)*(i-1))] ) ) 
        gamma <- c(gamma, unname( object$coefficients[(dim(object$x)[2]+1):(dim(object$x)[2]+ (m+2)) +(dim(object$x)[2]+m+2)*(i-1)] ))
      }
    }
  }
  
  ### type 
  
  if(type=="net") {
    
    ## survivalNET
    if ("dist" %in% names(object)) {
      
      ##pas de strate
      if (!("xlevels" %in% names(object))) {
        
        ##avec covariables
        if(dim(object$x)[2] != 0){
          net_cov <- function(x) { exp( exp(as.matrix(covariates)%*%beta)*(1-(1+(x/sigma)^nu)^(1/theta)) ) }
        }
        
        ##pas de covariables
        if(dim(object$x)[2] == 0){
          net_nocov <- function(x) { exp(1-(1+(x/sigma)^nu)^(1/theta) ) }
        }
        
      }
      
      ##avec strate
      if ("xlevels" %in% names(object)) {
        
        ##avec covariables
        if(dim(object$x)[2] != 0){
          
          net_strata_cov <- function(x, covariates, sigmas, nus, thetas) {
            n <- dim(covariates)[1]
            timpos <- dim(covariates)[2]
            rend <- data.frame()
            for (i in 1:n){
              timeval <- covariates[i,timpos]
              covi <- covariates[i,-timpos]
              sigmak <- as.vector(sigmas[timeval])
              nuk <- as.vector(nus[timeval])
              thetak <- as.vector(thetas[timeval])
              sur <-exp( exp(as.numeric(as.matrix(covi)%*%beta))*(1-(1+(x/sigmak)^nuk)^(1/thetak)) )
              rend <- rbind(rend, sur)
            }
            rend
          }
          
          
        }
        
        ##pas de covariables 
        if(dim(object$x)[2] == 0){
          
          net_strata_nocov <- function(x, covariates, sigmas, nus, thetas) {
            n <- dim(covariates)[1]
            timpos <- dim(covariates)[2]
            rend <- data.frame()
            for (i in 1:n){
              timeval <- covariates[i,timpos]
              sigmak <- sigmas[timeval]
              nuk <- nus[timeval]
              thetak <- thetas[timeval]
              sur <-exp( (1-(1+(x/sigmak)^nuk)^(1/thetak)) )
              rend <- rbind(rend, sur)
            }
            rend
          }  
          
        }
      }
    }
    
    ## survivalFLEXNET
    if("m" %in% names(object)){
      
      ##pas de strate  
      if (!("xlevels" %in% names(object))) {
        
        #valeurs de la spline
        splnvalues <- splinecube(newtimes, gamma, m, mpos )$spln
        ##avec covariables
        if(dim(object$x)[2] != 0){
          flex_net_cov <- function(x){ exp( -1*exp(as.matrix(covariates)%*%beta)*exp(x) )}
        }
        
        ##pas de covariables
        if(dim(object$x)[2] == 0){
          
          flex_net_nocov <- function(x){ exp( -1*exp(x) ) }
        }
        
      }
      
      ##avec strate 
      if ("xlevels" %in% names(object)) {
        
        ##avec covariables
        if(dim(object$x)[2] != 0){
          
          flex_net_strata_cov <- function(x, linpred){ exp(-1*exp(linpred)*exp(x)) }
        }
        
        ##pas de covariables 
        if(dim(object$x)[2] == 0){
          
          flex_net_strata_nocov <-  function(x){ exp( -1*exp(x) ) }
          
        }
      }  
      
    }
    
  }
  
  ###predictions 
  
  if (type == "net"){
    
    ## survivalNET
    if ("dist" %in% names(object)) {
      
      ##sans strate
      if (!("xlevels" %in% names(object))){
        
        ##avec covariables
        if(dim(object$x)[2] != 0){
          
          predictions <- matrix(sapply(newtimes, FUN = "net_cov"),nrow = n)
        }
        
        ##sans covariables
        if(dim(object$x)[2] == 0){
          
          predictions <- data.frame(matrix(rep(net_nocov(newtimes), each = n), nrow = n, byrow = FALSE))}
      }
      
      ##avec strate
      if ("xlevels" %in% names(object)){
        
        correstab <- object$correstab
        timevar <- as.character(unlist(covariates[names(object$xlevels)]))
        timevarnum <- as.numeric(correstab[timevar]) 
        covariates[,dim(covariates)[2]] <- timevarnum
        covariates <- data.frame(lapply(covariates, as.numeric))
        
        ##avec covariables
        if(dim(object$x)[2] != 0){
          
          predictions <- net_strata_cov(x = newtimes, covariates = covariates,
                                        sigmas = sigmas, nus = nus, thetas = thetas)
        }
        
        ##sans covariables
        if(dim(object$x)[2] == 0){
          
          predictions <- net_strata_nocov(x = newtimes, covariates = covariates,
                                          sigmas = sigmas, nus = nus, thetas = thetas)
        }
      }
      
    }
    
    ## survivalFLEXNET
    if ("m" %in% names(object)) {
      
      ##sans strate
      if (!("xlevels" %in% names(object))){
        
        ##avec covariables
        if(dim(object$x)[2] != 0){
          
          predictions <- matrix(sapply(splnvalues, FUN = "flex_net_cov"),nrow = n)
        }
        
        ##sans covariables 
        if(dim(object$x)[2] == 0){
          
          predictions <- data.frame(matrix(rep(sapply(splnvalues, FUN = "flex_net_nocov")
                                               , each = n), nrow = n, byrow = FALSE))}
      }
      
      ##avec strate
      if ("xlevels" %in% names(object)){
        
        ##avec covariables
        if(dim(object$x)[2] != 0){
          
          correstab <- object$correstab
          timevar <- as.character(unlist(covariates[names(object$xlevels)]))
          timevarnum <- as.numeric(correstab[timevar]) 
          K = sort(unique(timevarnum))
          covariates[,dim(covariates)[2]] <- timevarnum
          covariates <- data.frame(lapply(covariates, as.numeric))
          
          betas <-  matrix(beta, ncol = length(object$xlevels[[1]]))
          gammas <- matrix(gamma, ncol = length(object$xlevels[[1]]))
          
          n <- dim(covariates)[1]
          predictions <- as.data.frame(matrix(NA, nrow = n, ncol = length(newtimes)))
          
          for(k in K){
            
            gammak <- gammas[,k]
            betak <- betas[,k]
            mposk <- mpos[,k]
            idx <- which(timevarnum == k)    
            cov_k <- covariates[idx, -ncol(covariates), drop = FALSE]
            
            linpreds <- as.numeric(as.matrix(cov_k) %*% betak)
            splnvalues <- splinecube(newtimes, gammak, m, mposk )$spln
            
            surv_k <- sapply(splnvalues, FUN =  "flex_net_strata_cov", linpreds)
            
            predictions[idx, ] <- surv_k
          }
        }
        
        ##sans covariables
        if(dim(object$x)[2] == 0){
          
          correstab <- object$correstab
          timevar <- as.character(unlist(covariates[names(object$xlevels)]))
          timevarnum <- as.numeric(correstab[timevar]) 
          K = sort(unique(timevarnum))
          covariates[,dim(covariates)[2]] <- timevarnum
          covariates <- data.frame(lapply(covariates, as.numeric))
          
          gammas <- matrix(gamma, ncol = length(object$xlevels[[1]]))
          
          n <- dim(covariates)[1]
          predictions <- as.data.frame(matrix(NA, nrow = n, ncol = length(newtimes)))
          
          for(k in K){
            
            gammak <- gammas[,k]
            mposk <- mpos[,k]
            idx <- which(timevarnum == k)    
            
            splnvalues <- splinecube(newtimes, gammak, m, mposk )$spln
            
            surv_k <- sapply(splnvalues, FUN ="flex_net_strata_nocov")
            
            predictions[idx, ] <- as.data.frame(t(matrix(rep(surv_k, each = length(idx)), nrow = length(idx), byrow = TRUE) )) 
          }
        }
      } 
    }
    predictions <- unname(cbind(rep(1, dim(predictions)[1]), predictions))
  }
  
  if(type=="lp") { ##enlever la covar strate de covariates
    if ("xlevels" %in% names(object)){ 
      covariates <- covariates[, -which(names(covariates) == names(object$xlevels))]
    }
    
    predictions <-  covariates %*% beta 
  }
  
  newtimes <- c(0, newtimes)
  
  predictions <- as.data.frame(predictions)
  names(predictions) <- newtimes
  
  if(!is.null(newtimesSave)){
    predictions <- predictions[,c(as.character(newtimesSave))]
    
    newtimes <- newtimesSave
  }
  
  return(list(times=newtimes, predictions=predictions))
  
}