
predict.survivalNET <- function(object, type="relative", newdata=NULL, newtimes=NULL,
                                ratetable = NULL, ...){
    
  if(is.null(newtimes))  { newtimes <- 1:max(object$y[,1]) }
  
  if(0 %in% newtimes){
      newtimes <- sort(newtimes[-(newtimes == 0)])
    }
    
  if(!is.null(newdata))
    {
      if(!is.data.frame(newdata)) stop("Argument 'newdata' must be a data frame")
    
      covnames <- names(as.data.frame(object$x))
      if ("xlevels" %in% names(object)) {
        covnames <- c(covnames, names(object$xlevels)) 
      }
    
      indic <- covnames %in% names(newdata) 
      if( sum(!indic) > 0 ) stop("Missing predictor in the data frame")
      covariates <- as.matrix(newdata[,covnames])
      
      if(type=="overall"){
        if(is.null(ratetable))stop("For overall survival, the 'ratetable' argument is necessary.")
        all_terms = attr(terms(object$formula), "term.labels")
        ratetable_terms <- grep("ratetable\\(", all_terms, value = TRUE)
        extract_vars <- function(term) {
          var_string <- sub("^[^\\(]+\\((.*)\\)$", "\\1", term)
          vars <- trimws(unlist(strsplit(var_string, ",")))
          return(vars)
        }
        
        assign_ratetable_vars <- function(vars) {
          age <- year <- sex <- NULL
          for (var in vars) {
            if (grepl("age = ", var)) {
              age <- sub("age = ", "", var)
            } else if (grepl("year = ", var)) {
              year <- sub("year = ", "", var)
            } else if (grepl("sex = ", var)) {
              sex <- sub("sex = ", "", var)
            }
          }
          unnamed_vars <- setdiff(vars, c(age, sex, year))
          if (length(unnamed_vars) > 0) {
            if (is.null(age) && length(unnamed_vars) >= 1) age <- unnamed_vars[1]
            if (is.null(year) && length(unnamed_vars) >= 2) year <- unnamed_vars[2]
            if (is.null(sex) && length(unnamed_vars) >= 3) sex <- unnamed_vars[3]
          }
          return(list(age = age, sex = sex, year = year))
        }
        
        extract_vars <- function(term) {
          var_string <- sub("^[^\\(]+\\((.*)\\)$", "\\1", term)
          vars <- trimws(unlist(strsplit(var_string, ",")))
          return(vars)
        }
        
        ratetable_vars <- assign_ratetable_vars(unlist(lapply(ratetable_terms, extract_vars)))
        
        age = ratetable_vars$age
        year = ratetable_vars$year
        sex = ratetable_vars$sex
        indic <- c(age,year,sex, covnames)  %in% names(newdata) 
        if( sum(!indic) > 0 ) stop("Missing predictor in the data frame.
                                   For overall suvival, newdata also needs
                                   'age', 'sex' and 'year' for the ratetable")
      }
       }
    
  if(is.null(newdata))  { 
    
    covnames <- names(as.data.frame(object$x))
    covariates <- object$x 
    if ("xlevels" %in% names(object)) {
      covnames <- c(covnames, names(object$xlevels)) 
      covariates <- cbind(covariates, object$levelsval)
      colnames(covariates)[ncol(covariates)] <- names(object$xlevels)
      covariates <- as.data.frame(covariates)
      }
      
        if(type=="overall"){
          if(is.null(ratetable))stop("For overall survival, the 'ratetable' argument is necessary.")
          age = object$ays$age
          year = object$ays$year
          sex = object$ays$sex
        }
    
    }
    
    
  ### regression coefficients
    
  if ("dist" %in% names(object)) {
     
        if(object$dist=="genweibull")  {
              beta <- object$coefficients[1:(dim(object$x)[2])]
              sigma <- exp(object$coefficients[(dim(object$x)[2])+1])
              nu <- exp(object$coefficients[(dim(object$x)[2])+2])
              theta <- exp(object$coefficients[(dim(object$x)[2])+3])
              }
            
        if(object$dist=="weibull")  {
            beta <- object$coefficients[1:(dim(object$x)[2])]
            sigma <- exp(object$coefficients[(dim(object$x)[2])+1])
            nu <- exp(object$coefficients[(dim(object$x)[2])+2])
            theta <- 1
              }
            
        if(object$dist=="exponential")  {
            beta <- object$coefficients[1:(dim(object$x)[2])]
            sigma <- exp(object$coefficients[(dim(object$x)[2])+1])
            nu <- 1
            theta <- 1
              }
      
        if ("xlevels" %in% names(object)) {   
        
        if(object$dist=="genweibull")  {
          beta <- object$coefficients[1:(dim(object$x)[2])]
          sigmas <- exp(object$coefficients[grep("^log sigma_", 
                                names(object$coefficients))])
          nus <- exp(object$coefficients[grep("^log nu_", 
                                names(object$coefficients))])
          thetas <- exp(object$coefficients[grep("^log theta_", 
                                names(object$coefficients))])
        }
        if(object$dist=="weibull")  {
          beta <- object$coefficients[1:(dim(object$x)[2])]
          sigmas <- exp(object$coefficients[grep("^log sigma_", 
                                names(object$coefficients))])
          nus <- exp(object$coefficients[grep("^log nu_", 
                                names(object$coefficients))])
          thetas <- rep(1,length(sigmas))
        }
        if(object$dist=="exponential")  {
          beta <- object$coefficients[1:(dim(object$x)[2])]
          sigmas <- exp(object$coefficients[grep("^log sigma_", 
                             names(object$coefficients))])
          nus <- rep(1,length(sigmas))
          thetas <- rep(1,length(sigmas))
        }
       }
    }
  
  if("m" %in% names(object)){
      
       beta <- object$coefficients[1:(dim(object$x)[2])]
       gamma <- object$coefficients[(dim(object$x)[2]+1):(dim(object$t.table)[1])]
       m = object$m
       mpos = object$mpos
  }
  
  ### type 
    
  if(type=="relative") {
    
    if ("dist" %in% names(object)) {
      
      if (!("xlevels" %in% names(object))) {
      
      if(dim(object$x)[2] != 0){
      fun <- function(x) { exp( exp(covariates%*%beta)*(1-(1+(x/sigma)^nu)^(1/theta)) ) }
      }
      
      if(dim(object$x)[2] == 0){
        fun <- function(x) { exp(1-(1+(x/sigma)^nu)^(1/theta) ) }
      }
      
      }
      
      if ("xlevels" %in% names(object)) {
         
        if(dim(object$x)[2] != 0){
           
          # fun <- function(x, covariates, sigmas, nus, thetas, timecov, K) {
          #   rend <- data.frame()
          #   for (k in K){
          #     covariatesk <- as.matrix( covariates[covariates[[timecov]] == k, -dim(covariates)[2]]) 
          #     sigmak <- sigmas[k]
          #     nuk <- nus[k]
          #     thetak <- thetas[k]
          #     sur <-exp( exp(covariatesk%*%beta)*(1-(1+(x/sigmak)^nuk)^(1/thetak)) )
          #   rend <- rbind(rend, sur)
          #   }
          #   rend[order(rownames(rend)),]
          # }
          
          fun <- function(x, covariates, sigmas, nus, thetas, timecov, K) {
            n <- dim(covariates)[1]
            timpos <- dim(covariates)[2]
            rend <- data.frame()
            for (i in 1:n){
              timeval <- covariates[i,timpos]
              covi <- covariates[i,-timpos]
              sigmak <- sigmas[timeval]
              nuk <- nus[timeval]
              thetak <- thetas[timeval]
              sur <-exp( exp(covi%*%beta)*(1-(1+(x/sigmak)^nuk)^(1/thetak)) )
              rend <- rbind(rend, sur)
            }
            rend
          }
          
          
        }
        
        if(dim(object$x)[2] == 0){
          
          fun <- function(x, covariates, sigmas, nus, thetas, timecov, K) {
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
    
    if("m" %in% names(object)){
      
      if (!("xlevels" %in% names(object))) {
      
        splnvalues <- splinecube(newtimes, gamma, m, mpos)$spln
      
      if(dim(object$x)[2] != 0){
        fun <- function(x){ exp( -1*exp(as.matrix(covariates)%*%beta)*exp(x))}
      }
      
      if(dim(object$x)[2] == 0){
        fun <- function(x){ exp( -1*exp(x) ) }
      }
      
      }
      
      if ("xlevels" %in% names(object)) {
      
        if(dim(object$x)[2] != 0){
          
          gammas <- matrix(object$coefficients[-(1:dim(object$x)[2])], ncol = length(object$xlevels[[1]]))
          
          fun <- function(x, covariates, gammas, timecov, K=NULL) {
            n <- dim(covariates)[1]
            timpos <- dim(covariates)[2]
            rend <- data.frame()
            for (i in 1:n){
              timeval <- covariates[i,timpos]
              gammai <- gammas[,timeval]
              covariatesi <- covariates[i,-timpos]
              splnvalues <- splinecube(x, gammai, m, mpos)$spln
              sur <-exp(-1*exp(as.matrix(covariatesi)%*%beta)*exp( splnvalues ))
              rend <- rbind(rend, sur)
            }
            rend
          }
        }
        
        if(dim(object$x)[2] == 0){
          
          gammas <- matrix(object$coefficients, ncol = length(object$xlevels[[1]]))
          
          fun <- function(x, covariates, gammas, timecov, K=NULL) {
            n <- dim(covariates)[1]
            timpos <- dim(covariates)[2]
            rend <- data.frame()
            for (i in 1:n){
              timeval <- covariates[i,timpos]
              gammai <- gammas[,timeval]
              splnvalues <- splinecube(x, gammai, m, mpos)$spln
              sur <-exp(-1*exp( splnvalues ))
              rend <- rbind(rend, sur)
            }
            rend
          }
        }
   }  
    
  }
  }
  
  if(type=="lp") { ##enlever la covar strate de covariates pour ca 
    fun <- function(x) { covariates %*% beta }
    }
    
  if(type=="overall"){
    if ("dist" %in% names(object)) {
      fun <- function(x) {
        sapply(x, function(t) {
          exp(exp(covariates %*% beta) * (1 - (1 + (t / sigma)^nu)^(1 / theta))) * 
            exp(-1 * sapply(1:dim(covariates)[1], 
                            FUN = function(i) {
                              expectedcumhaz(ratetable, age[i], year[i], sex[i], t)
                            }))
        })
      }
    }
    if("m" %in% names(object)){
      splnvalues <- splinecube(newtimes, gamma, m, mpos)$spln
      fun <- function(x) {
        sapply(x, function(t) {
          exp(- exp(covariates %*% beta + x) ) * 
            exp(-1 * sapply(1:dim(covariates)[1], 
                            FUN = function(i) {
                              expectedcumhaz(ratetable, age[i], year[i], sex[i], t)
                            }))
        })
      }
    } 
  }
      
      
  ###predictions 
  
  if ("dist" %in% names(object)) {
      
      if (!("xlevels" %in% names(object))){
        
        if(dim(object$x)[2] != 0){
             predictions <- sapply(newtimes, FUN = "fun")}
        
        if(dim(object$x)[2] == 0){
             n <- dim(object$y)[1]
             predictions <- data.frame(matrix(rep(fun(newtimes), each = n), nrow = n, byrow = FALSE))}
      }
   
      if ("xlevels" %in% names(object)){
       
         correstab <- object$correstab
         timevar <- as.integer(unlist(covariates[names(object$xlevels)]))
         timevarnum <- as.numeric(correstab[as.character(timevar)]) 
         covariates[,dim(covariates)[2]] <- timevarnum
         timecov <- names(object$xlevels)
           
         K = sort(unique(timevarnum))
  
         predictions <- fun(x = newtimes, covariates = covariates,
                            sigmas = sigmas, nus = nus, thetas = thetas,
                            timecov = timecov, K = K)
     }
    
      }
  
  if ("m" %in% names(object)) {
      
      if (!("xlevels" %in% names(object))){
        
          if(dim(object$x)[2] != 0){
              predictions <- sapply(splnvalues, FUN = "fun")}
          if(dim(object$x)[2] == 0){
              n <- dim(object$y)[1]
              predictions <- data.frame(matrix(rep(sapply(splnvalues, FUN = "fun")
                                        , each = n), nrow = n, byrow = FALSE))}
         }
    
      if ("xlevels" %in% names(object)){
          
        correstab <- object$correstab
        timevar <- as.integer(unlist(covariates[names(object$xlevels)]))
        timevarnum <- as.numeric(correstab[as.character(timevar)]) 
        covariates[,dim(covariates)[2]] <- timevarnum
        timecov <- names(object$xlevels)
        
        K = sort(unique(timevarnum))
        
          if(dim(object$x)[2] != 0){
              n <- dim(object$y)[1]
              predictions <- fun(newtimes, covariates = 
                                 covariates, gammas = gammas, timecov = 
                                 timecov, K= K)}
          
          if(dim(object$x)[2] == 0){
              n <- dim(object$y)[1]
              predictions <- fun(newtimes, covariates = 
                                  covariates, gammas = gammas, timecov = 
                                  timecov, K= K)}
      }
        } 
  
   predictions <- unname(cbind(rep(1, dim(predictions)[1]), predictions))
  
  if(type == "lp"){
    predictions <- predictions[,2]
  }
  
  newtimes <- c(0, newtimes)
  
  return(list(times=newtimes, predictions=predictions))
  }
  
  
  
  
  
