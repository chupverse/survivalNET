
predict.survivalNET <- function(object, type="net", newdata=NULL, newtimes=NULL,
                                ratetable = NULL, method = NULL, ...){
    
  if(!(type %in% c("net","lp","overall")))  stop("Argument 
                  'type' must be 'net', 'lp' or 'overall' ")
  if(type == "overall" && "m" %in% names(object))stop("The 'overall' survival prediction
  for survivalFLEXNET is still under development. Please
                            use 'type = 'relative' or 'lp'. ")
  
  if(is.null(newtimes))  { 
    newtimes <- 1:max(object$y[,1])
                           
    newtimesSave <- NULL
    
  }else{
    
    newtimesSave <- newtimes
        
    newtimes <- sort(c((1:max(object$y[,1])),unique(newtimes)))}
  
  if(0 %in% newtimes){
      newtimes <- sort(newtimes[-(newtimes == 0)])
    }
    
  if(!is.null(newdata))
    { 
    n <- dim(newdata)[1]
      if(!is.data.frame(newdata)) stop("Argument 'newdata' must be a data frame")
    
      covnames <- names(as.data.frame(object$x))
      if ("xlevels" %in% names(object)) {
        covnames <- c(covnames, names(object$xlevels)) 
      }
    
      indic <- covnames %in% names(newdata) 
      if( sum(!indic) > 0 ) stop("Missing predictor in the data frame")
      covariates <- as.data.frame(as.matrix(newdata[,covnames]))
      names(covariates) <- covnames
      
      if(type=="overall"){
        if(is.null(ratetable))stop("For overall survival, the 'ratetable' argument is necessary.")
        if(!is.null(method) && !(method %in% c("exact", "trapezoidal", "table")))stop("For overall 
                          survival, method can only take the values 'exact', 'trapzoidal' or 
                          'table'.")
        if(is.null(method)){
          method = "exact"}else{
          method = method
        }
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
    
    n <- dim(object$y)[1]
    
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
          if(!is.null(method) && !(method %in% c("exact", "trapezoidal", "table")))stop("For overall 
                          survival, method can only take the values 'exact', 'trapzoidal' or 
                          'table'.")
          if(is.null(method)){
            method = "exact"}else{
              method = method
            }
          age = object$ays$age
          year = object$ays$year
          sex = object$ays$sex
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
      
       beta <- unname( object$coefficients[1:(dim(object$x)[2])] )
       gamma <- unname( object$coefficients[(dim(object$x)[2]+1):(dim(object$t.table)[1])] ) 
       m <- object$m
       mpos <- object$mpos
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
                
                net_strata_cov <- function(x, covariates, sigmas, nus, thetas, timecov, K) {
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
            
            ##pas de covariables 
              if(dim(object$x)[2] == 0){
                
                net_strata_nocov <- function(x, covariates, sigmas, nus, thetas, timecov, K) {
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
                
                gammas <- matrix(object$coefficients[-(1:dim(object$x)[2])], ncol = length(object$xlevels[[1]]))
                
                flex_net_strata_cov <- function(x, covariates, gammas, timecov, K=NULL) {
                  n <- dim(covariates)[1]
                  timpos <- dim(covariates)[2]
                  rend <- data.frame()
                  for (i in 1:n){
                    timeval <- covariates[i,timpos]
                    gammai <- gammas[,timeval]
                    covariatesi <- covariates[i,-timpos]
                    splnvalues <- splinecube(x, gammai, m, mpos )$spln
                    sur <-exp(-1*exp(as.vector(as.numeric(covariatesi)%*%beta))*exp( splnvalues ))
                    rend <- rbind(rend, sur)
                  }
                  rend
                }
              }
             
            ##pas de covariables 
              if(dim(object$x)[2] == 0){
                
                gammas <- matrix(object$coefficients, ncol = length(object$xlevels[[1]]))
                
                flex_net_strata_nocov <- function(x, covariates, gammas, timecov, K=NULL) {
                  n <- dim(covariates)[1]
                  timpos <- dim(covariates)[2]
                  rend <- data.frame()
                  for (i in 1:n){
                    timeval <- covariates[i,timpos]
                    gammai <- gammas[,timeval]
                    splnvalues <- splinecube(x, gammai, m, mpos )$spln
                    sur <-exp(-1*exp( splnvalues ))
                    rend <- rbind(rend, sur)
                  }
                  rend
                }
              }
       }  
      
    }
  
    }
  
  if(type=="overall"){
    
    ## survivalNET
      if ("dist" %in% names(object)) {
        
      ## pas de strate
        if (!("xlevels" %in% names(object))){
          
          ##avec covariables
              if(dim(object$x)[2] != 0){    
            
                net_ov_cov <- function(x) {
              sapply(x, function(t) {
                exp(exp(covariates %*% beta) * (1 - (1 + (t / sigma)^nu)^(1 / theta))) * 
                  exp(-1 * sapply(1:dim(covariates)[1], 
                                  FUN = function(i) {
                                    expectedcumhaz(ratetable, object$ays$age[i], 
                                                   object$ays$year[i], object$ays$sex[i], t,
                                                   method = method)
                                                    }))
                                    })
                                }
            }
           
          ##pas de covariables 
              if(dim(object$x)[2] == 0){
              
                net_ov_nocov <- function(x) {
                sapply(x, function(t) {
                  exp((1 - (1 + (t / sigma)^nu)^(1 / theta))) * 
                    exp(-1 * sapply(1:dim(covariates)[1], 
                                    FUN = function(i) {
                                      expectedcumhaz(ratetable, object$ays$age[i], 
                                                     object$ays$year[i], object$ays$sex[i], t, 
                                                     method = method)
                                    }))
                })
              }
            }
        
        }
        
      ## avec strate 
        if ("xlevels" %in% names(object)){
          
          ##avec covariables
              if(dim(object$x)[2] != 0){
                net_ov_strata_cov <- function(x, covariates, sigmas, nus, thetas, timecov, K) {
                  n <- dim(covariates)[1]
                  timpos <- dim(covariates)[2]
                  rend <- data.frame()
                  for (i in 1:n){
                    timeval <- covariates[i,timpos]
                    covi <- covariates[i,-timpos]
                    sigmak <- sigmas[timeval]
                    nuk <- nus[timeval]
                    thetak <- thetas[timeval]
                    expected_values <- sapply(x, 
                              FUN = function(t){expectedcumhaz(ratetable, object$ays$age[i], 
                                                            object$ays$year[i], object$ays$sex[i], 
                                                            t, method = method)})
                    sur <- exp( exp(as.vector(as.numeric(covi)%*%beta))*(1-(1+(x/sigmak)^nuk)^(1/thetak)))* 
                      exp(-1 * expected_values)
                    rend <- rbind(rend, sur)
                  }
                  rend
                }
              }
            
          ##pas de covariables
              if(dim(object$x)[2] == 0){
               
                 net_ov_strata_nocov <- function(x, covariates, sigmas, nus, thetas, timecov, K) {
                  n <- dim(covariates)[1]
                  timpos <- dim(covariates)[2]
                  rend <- data.frame()
                  for (i in 1:n){
                    timeval <- covariates[i,timpos]
                    sigmak <- sigmas[timeval]
                    nuk <- nus[timeval]
                    thetak <- thetas[timeval]
                    expected_values <- sapply(x, 
                                    FUN = function(t){expectedcumhaz(ratetable, object$ays$age[i], 
                                                            object$ays$year[i], object$ays$sex[i], 
                                         t, method = method)})
                    sur <-exp( (1-(1+(x/sigmak)^nuk)^(1/thetak)))* 
                      exp(-1 * expected_values)
                    rend <- rbind(rend, sur)
                  }
                  rend
                 } ##fin fonction
                 
              }
        }
      }
    
    ## survivalFLEXNET
      if("m" %in% names(object)){
  
        #valeurs de la spline
        splnvalues <- splinecube(newtimes, gamma, m, mpos )$spln
       
        ## pas de strate    
          if (!("xlevels" %in% names(object))){
              
            ##avec covariables
              if(dim(object$x)[2] != 0){
              
            flex_ov_cov <- function(x) {
              sapply(x, function(t, newtimes, splnvalues) {
                exp(- exp(splnvalues[which(newtimes == t)] + as.matrix(covariates) %*% beta)) *
                  exp(-1 * sapply(1:dim(covariates)[1], 
                        FUN = function(i){
                        expectedcumhaz(ratetable, object$ays$age[i], 
                                       object$ays$year[i], object$ays$sex[i], t)}))
              }, newtimes = newtimes, splnvalues = splnvalues)
             }
            }
           
            ##pas de covariables  
              if(dim(object$x)[2] == 0){
                
                flex_ov_nocov <- function(x) {
                  sapply(x, function(t) {
                    exp(- exp( x ) ) * 
                      exp(-1 * sapply(1:dim(covariates)[1], 
                                      FUN = function(i){
                                        expectedcumhaz(ratetable, object$ays$age[i], 
                                                       object$ays$year[i], object$ays$sex[i], t)}))
                  })
                }
              }
              
          }
         
        ## avec strate           
          if ("xlevels" %in% names(object)){
                
            ##avec covariables
              if(dim(object$x)[2] != 0){
                  gammas <- matrix(object$coefficients[-(1:dim(object$x)[2])],
                                   ncol = length(object$xlevels[[1]]))
                  
                  flex_ov_strata_cov <- function(x, covariates, timecov, K) {
                    n <- dim(covariates)[1]
                    timpos <- dim(covariates)[2]
                    rend <- data.frame()
                    for (i in 1:n){
                      timeval <- covariates[i,timpos]
                      covi <- covariates[i,-timpos]
                      gammai <- gammas[,timeval]
    
                      expected_values <- sapply(x, 
                                   FUN = function(t){
                                            expectedcumhaz(ratetable, object$ays$age[i], 
                                            object$ays$year[i], object$ays$sex[i], 
                                            t, method = method)})
                      sur <- exp( exp(as.vector(covi%*%beta))   ) * 
                        exp(-1 * expected_values)
                      rend <- rbind(rend, sur)
                    }
                    rend
                  }
                }
            
            ##pas de covariables
              if(dim(object$x)[2] == 0){
                  flex_ov_strata_nocov <- function(x, covariates, timecov, K){} 
                  }
        } 
    }
  
    } 
  
  if(type=="lp") { ##enlever la covar strate de covariates pour ca 
    if ("xlevels" %in% names(object)){ 
      covariates1 <- covariates[, -which(names(covariates) == names(object$xlevels))]
    }
    
      fun <- function(x) { covariates1 %*% beta }
  }
  
  ###predictions 
  
  if (type == "net"){
  
    ## survivalNET
      if ("dist" %in% names(object)) {
          
        ##sans strate
          if (!("xlevels" %in% names(object))){
            
            ##avec covariables
              if(dim(object$x)[2] != 0){
                
                   predictions <- sapply(newtimes, FUN = "net_cov")}
            
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
                timecov <- names(object$xlevels)
                
                K = sort(unique(timevarnum))
            
             ##avec covariables
               if(dim(object$x)[2] != 0){
               
                 predictions <- net_strata_cov(x = newtimes, covariates = covariates,
                                  sigmas = sigmas, nus = nus, thetas = thetas,
                                  timecov = timecov, K = K)
               }
              
             ##sans covariables
               if(dim(object$x)[2] == 0){
                 
                 predictions <- net_strata_nocov(x = newtimes, covariates = covariates,
                                     sigmas = sigmas, nus = nus, thetas = thetas,
                                     timecov = timecov, K = K)
                 }
          }
        
      }
      
    ## survivalFLEXNET
      if ("m" %in% names(object)) {
          
        ##sans strate
          if (!("xlevels" %in% names(object))){
            
             ##avec covariables
                if(dim(object$x)[2] != 0){
                    
                  predictions <- sapply(splnvalues, FUN = "flex_net_cov")}
             
             ##sans covariables 
                if(dim(object$x)[2] == 0){
                    
                    predictions <- data.frame(matrix(rep(sapply(splnvalues, FUN = "flex_net_nocov")
                                              , each = n), nrow = n, byrow = FALSE))}
             }
        
        ##avec strate
          if ("xlevels" %in% names(object)){
              
            correstab <- object$correstab
            timevar <- as.character(unlist(covariates[names(object$xlevels)]))
            timevarnum <- as.numeric(correstab[timevar]) 
            covariates[,dim(covariates)[2]] <- timevarnum
            timecov <- names(object$xlevels)
            
            K = sort(unique(timevarnum))
              
              ##avec covariables
                if(dim(object$x)[2] != 0){
                    
                    predictions <- flex_net_strata_cov(newtimes, covariates = 
                                       covariates, gammas = gammas, timecov = 
                                       timecov, K= K)
                    }
              
              ##sans covariables
                if(dim(object$x)[2] == 0){
                    
                    predictions <- flex_net_strata_nocov(newtimes, covariates = 
                                        covariates, gammas = gammas, timecov = 
                                        timecov, K= K)}
          }
  } 
    
  predictions <- unname(cbind(rep(1, dim(predictions)[1]), predictions))
    
  }
  
  if (type == "overall"){
    
    ## survivalNET
      if ("dist" %in% names(object)) {
        
        ##sans strate
          if (!("xlevels" %in% names(object))){
            
            ##avec covariables
                if(dim(object$x)[2] != 0){
                  
                  predictions <- sapply(newtimes, FUN = "net_ov_cov")}
              
            ##sans covariables
                if(dim(object$x)[2] == 0){
                  
                  predictions <- net_ov_nocov(newtimes)}
          }
        
        ##avec strate
          if ("xlevels" %in% names(object)){
            
            correstab <- object$correstab
            timevar <- as.character(unlist(covariates[names(object$xlevels)]))
            timevarnum <- as.numeric(correstab[timevar]) 
            covariates[,dim(covariates)[2]] <- timevarnum
            timecov <- names(object$xlevels)
            
            K = sort(unique(timevarnum))
              
              ##avec covariables
                if(dim(object$x)[2] != 0){
                  
                  predictions <- net_ov_strata_cov(x = newtimes, covariates = covariates,
                                    sigmas = sigmas, nus = nus, thetas = thetas,
                                    timecov = timecov, K = K)
                }
              
              ##sans covariables
                if(dim(object$x)[2] == 0){
                  
                  predictions <- net_ov_strata_nocov(x = newtimes, covariates = covariates,
                                       sigmas = sigmas, nus = nus,
                                       thetas = thetas, timecov = timecov, K = K)
                }
          }
        
      }
      
    ## survivalFLEXNET
      if ("m" %in% names(object)) {
        
        ##sans strate
          if (!("xlevels" %in% names(object))){
            
             ##avec covariables
                if(dim(object$x)[2] != 0){
                  predictions <- flex_ov_cov(newtimes)
                }
             ##sans covariables
                if(dim(object$x)[2] == 0){

                  predictions <- sapply(splnvalues, FUN = "flex_ov_nocov")
                                                       }
          }
          
        ##avec strate
          if ("xlevels" %in% names(object)){
            
            correstab <- object$correstab
            timevar <- as.character(unlist(covariates[names(object$xlevels)]))
            timevarnum <- as.numeric(correstab[timevar]) 
            covariates[,dim(covariates)[2]] <- timevarnum
            timecov <- names(object$xlevels)
            
            K = sort(unique(timevarnum))
                
              ##avec covariables  
                if(dim(object$x)[2] != 0){

                      predictions <- flex_ov_strata_cov(newtimes, covariates = 
                                       covariates, gammas = gammas, timecov = 
                                       timecov, K= K)}
                
              ##sans covariables
                if(dim(object$x)[2] == 0){
                  
                  predictions <- flex_ov_strata_nocov(newtimes, covariates = 
                                       covariates, gammas = gammas, timecov = 
                                       timecov, K= K)}
          }
      } 
      
    predictions <- unname(cbind(rep(1, dim(predictions)[1]), predictions))
    
  }
   
  if(type == "lp"){
    predictions <- covariates %*% object$coefficients[covnames]
  }
  
  newtimes <- c(0, newtimes)
  
  predictions <- as.data.frame(predictions)
  names(predictions) <- newtimes
  
  if(!is.null(newtimesSave)){
    predictions <- predictions[,c(as.character(newtimesSave)) ]
    
    newtimes <- newtimesSave
  }
  
  return(list(times=newtimes, predictions=predictions))
  
  }
  
  
  
  
  
