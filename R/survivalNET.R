
survivalNET <- function(formula, data, ratetable, dist="weibull", weights=NULL)

{
  
  call = match.call()
  
  ####### check errors
  
  if (missing(formula)) stop("a formula argument is required")
  if (missing(data)) stop("a data argument is required")
  if (missing(ratetable)) stop("a ratetable argument is required")
  if (as.character(class(formula)) != "formula") stop("The first argument must
                   be a formula")
  if (as.character(class(data)) != "data.frame") stop("The second argument must
                   be a data frame")
  if (length(dim(ratetable))!=3) stop("The life table must have 3 dimensions:
                   age, year, sex")
  if (dim(ratetable)[3]!=2) stop("The life table must have 3 dimensions: age,
                   year, sex")
  if(!(dist %in% c("exponential","weibull","genweibull")))  stop("Argument 
                  'dist' must be 'exponential', 'weibull', or 'genweibull' ")
  
  ####### data management
  
  time <- data[,as.character(formula[[2]][2])] 
  event <- data[,as.character(formula[[2]][3])]
  all_terms <- attr(terms(formula), "term.labels")
  strata_terms <- grep("strata\\(", all_terms, value = TRUE)
  if(length(strata_terms)>1) stop("More than one 'strata' term found in  the formula. Only one variable at a time can be stratified")
  ratetable_terms <- grep("^ratetable\\(", all_terms, value = TRUE)
  if(length(ratetable_terms) == 0) stop("Error: The formula must contain a ratetable() term.")
  if(length(ratetable_terms)>1) stop("More than one 'ratetable' term found in  the formula.")
  CV <- setdiff(all_terms, c(strata_terms, ratetable_terms))
  if(length(CV) == 0){covnames = "1"} else{covnames <- CV}
  label<- NULL
  covs <- as.formula(paste("~", paste(covnames, collapse = " + ")))
  cova <- model.matrix(covs, data)[, -1, drop = FALSE]
  covnames <- colnames(cova)
  
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
  
  strata_var = unlist(lapply(strata_terms, extract_vars))
  
  if(!is.null(strata_var) && strata_var %in% covnames) stop("The stratified covariate also appears as a covariate in the formula.")
  
  if(is.null(strata_var)){
    timevar = strata_var
    xlevels = NULL
  }
  if(!is.null(strata_var)){
    timevar <- data[,strata_var]
    xlevels <- list(levels(as.factor(timevar)))
    names(xlevels) <- c(strata_var)
  }
  if(!is.null(timevar)){
    if(length(unique(timevar))>15) stop("The variable with a time-dependant effect has too many categories (>15)")
  }
  ratetable_vars <- assign_ratetable_vars(unlist(lapply(ratetable_terms, extract_vars)))
  age <- data[,ratetable_vars$age] 
  year <- data[,ratetable_vars$year]
  sex <- data[,ratetable_vars$sex] 
  if(is.null(weights)){weights = rep(1,dim(data)[1])}
  if(max(time)<30) stop("The event times must be expressed in days (max(", as.character(formula[[2]][2]),") is less than 30 days).")
  if(!is.numeric(age)) stop("'age' must be numeric")
  if(max(age)<30) stop("The ages must be expressed in days (max(", as.character(formula[[2]][2]),") is less than 30 days).")
  if(!is.character(sex)) stop("'sex' must be a character string")
  if(min(names(table(sex)) %in% c("female","male", NA))==0) stop("Argument 'sex' must be 'male' or 'female'")
  if(!is.date(year)) stop("The values for 'year' must be of the 'date' class")
  if(!is.null(weights)){
    if(!is.numeric(weights)) stop("Argument 'weights' must be a numeric vector")
    if(length(weights)!=dim(data)[1]) stop("Argument 'weights' must have the same length as the number of rows of the 'data' argument. (", dim(data)[1],")")}
  

  d <- cbind(time, event, cova, age, 1*(sex=="male"), year)
  na <- !is.na(apply(d, MARGIN=1, FUN = "sum"))
  
  time <- time[na]
  event <- event[na]
  cova <- as.matrix(cova[na,])
  age <- age[na]
  sex <- sex[na]
  year <- year[na]
  
  ###### Compute the expected mortality of individuals
  
  hP <- sapply(seq_along(time), function(i) {
    expectedhaz(ratetable = ratetable, age = age[i], sex = sex[i], year = year[i], time = time[i])
  })
  
  ###### log likelihood functions
  
  if (!is.null(covnames) & is.null(timevar)){
    
    
    logll1 <- function(sigma, nu, theta, beta, time, event, cova, hP, w){
      return(-1*sum( w*( event*log(hP+exp(cova%*%beta)*(
        (1/theta)*(1+(time/sigma)^nu)^((1/theta)-1)*(nu/sigma)*(time/sigma)^(nu-1) ) )
        + exp(cova%*%beta)*(1-(1+(time/sigma)^nu)^(1/theta)) ) ) )}
    
    if (dist=="genweibull"){
      
      init1 <- c(rep(0,dim(cova)[2]),rep(0,3))
      
      loglik1 <- function(par, time, event, cova, hP, w){
        dimC = dim(cova)[2]
        beta <- par[1:dimC]
        sigma <- exp(par[(dimC+1)]); nu <- exp(par[(dimC+2)]);  theta <- exp(par[(dimC+3)]);  
        return(logll1(sigma, nu, theta, beta, time, event, cova, hP, w)) }
      
    }
    
    if (dist=="weibull"){
      
      init1 <- c(rep(0,dim(cova)[2]),rep(0,2))
      
      loglik1 <- function(par, time, event, cova, hP, w){
        dimC = dim(cova)[2]
        beta <- par[1:dimC]
        sigma <- exp(par[(dimC+1)]); nu <- exp(par[(dimC+2)]);  theta <- 1 
        return(logll1(sigma, nu, theta, beta, time, event, cova, hP, w)) }
      
    }
    
    if (dist=="exponential"){
      
      init1 <- c(rep(0,dim(cova)[2]),rep(0,1))
      
      loglik1 <- function(par, time, event, cova, hP, w){
        dimC = dim(cova)[2]
        beta <- par[1:dimC]
        sigma <- exp(par[(dimC+1)]); nu <- 1;  theta <- 1
        return(logll1(sigma, nu, theta, beta, time, event, cova, hP, w)) }
    }
    
    suppressWarnings({
      logllmax1 <- optim(par = init1, fn = loglik1, time = time, event = event,
                         cova = cova, hP = hP, w = weights)
      
      indic <- 0
      while(indic <= 5){
        ll_val <- logllmax1$value
        logllmax1 <- optim(par = logllmax1$par, fn = loglik1, time = time, event = event,
                           cova = cova, hP = hP, w = weights)
        delta <- ll_val - logllmax1$value
        if(delta ==0) {indic = indic + 1}
      }
      
      logllmax1 <- optim(par = logllmax1$par, fn = loglik1, time = time, event = event,
                         cova = cova, hP = hP, , w = weights, hessian = TRUE)
    })  
    
  }
  
  if (!is.null(covnames) & !is.null(timevar)){
    
    timevarnames <- sort(unique(timevar))
    correstab <- setNames(seq_along(timevarnames), timevarnames)
    
    timevarnum <- as.numeric(correstab[as.character(timevar)])
    
    K = sort(unique(timevarnum))
    
    logll2 <- function(sigma, nu, theta, beta, time, event, cova, covatime, hP, w, K){
      
      value = 0
      for(k in K){
        sigmak <- sigma[k]
        nuk <- nu[k]
        thetak <- theta[k]
        timek <- time[covatime == k]
        eventk <- event[covatime == k]
        hPk <- hP[covatime == k]
        covak <- cova[covatime == k]
        wk <- w[covatime == k]
        value_strate <- -1*sum( wk*( eventk*log(hPk+exp(as.matrix(covak)%*%beta)*(
          (1/thetak)*(1+(timek/sigmak)^nuk)^((1/thetak)-1)*(nuk/sigmak)*(timek/sigmak)^(nuk-1) ) )
          + exp(as.matrix(covak)%*%beta)*(1-(1+(timek/sigmak)^nuk)^(1/thetak)) ) )
        value <- value +value_strate
      }
      return(value)
      
    }
    
    if (dist=="genweibull"){
      for (i in timevarnames){
        assign(paste0("param_names", i), c(paste0("log sigma_",i), paste0("log nu_",i), paste0("log theta_",i)))}
      
      init1 <- c(rep(0,dim(cova)[2]), rep(0,3*length(K))) 
      
      loglik2 <- function(par, time, event, cova, covatime, hP, w, K){
        dimC <- dim(cova)[2]
        beta <- par[1:dimC]
        sigma <- exp(par[(dimC+1):(dimC+length(K))])
        nu <- exp(par[(dimC+length(K)+1):(dimC+2*length(K))]) 
        theta <- exp(par[(dimC+2*length(K)+1):(dimC+3*length(K))])  
        return(logll2(sigma, nu, theta, beta, time, event, cova, covatime, hP, w, K)) }
      
      label <- covnames
      for (j in 1:3){
        for (i in timevarnames) {
          label <- c(label, get(paste0("param_names", i))[j] ) }
      }
    }
    
    if (dist=="weibull"){
      for (i in timevarnames){
        assign( paste0("param_names", i), c(paste0("log sigma_",i), paste0("log nu_",i)))}
      
      init1 <- c(rep(0,dim(cova)[2]), rep(0,2*length(K)))
      
      loglik2 <- function(par, time, event, cova, covatime, hP, w, K){
        dimC <- dim(cova)[2]
        beta <- par[1:dimC]
        sigma <- exp(par[(dimC+1):(dimC+length(K))])
        nu <- exp(par[(dimC+length(K)+1):(dimC+2*length(K))])
        theta <- rep(1,length(K)) 
        return(logll2(sigma, nu, theta, beta, time, event, cova, covatime, hP, w,K)) }
      
      label <- covnames
      for (j in 1:2){
        for (i in timevarnames) {
          label <- c(label, get(paste0("param_names", i))[j] ) }
      }
    }
    
    if (dist=="exponential"){
      for (i in timevarnames){
        assign(paste0("param_names", i), c(paste0("log sigma_",i)))}
      
      init1 <- c(rep(0,dim(cova)[2]),rep(0,length(K)))
      
      loglik2 <- function(par, time, event, cova, covatime, hP, w, K){
        dimC <- dim(cova)[2]
        beta <- par[1:dimC]
        sigma <- exp(par[(dimC+1):(dimC+length(K))])
        nu <- rep(1,length(K))
        theta <- rep(1,length(K))  
        return(logll2(sigma, nu, theta, beta, time, event, cova, covatime, hP, w, K)) }
      
      label <- covnames
      for (i in timevarnames) {
        label <- c(label, get(paste0("param_names", i)))}
      
    }
    
    suppressWarnings({
      logllmax1 <- optim(par = init1, fn = loglik2, time = time, event = event,
                         cova = cova, covatime = timevarnum, hP = hP, w = weights
                         , K= K)
      
      indic <- 0
      while(indic <= 5){
        ll_val <- logllmax1$value
        logllmax1 <- optim(par = logllmax1$par, fn = loglik2, time = time, 
                           event = event, cova = cova, covatime = timevarnum,
                           hP = hP, w = weights, K = K)
        delta <- ll_val - logllmax1$value
        if(delta ==0) {indic = indic + 1}
      }
      
      logllmax1 <- optim(par = logllmax1$par, fn = loglik2, time = time, 
                         event = event, cova = cova, covatime = timevarnum,
                         hP = hP, w = weights, K= K,
                         hessian = TRUE)
    })
  }
  
  if (is.null(covnames) & !is.null(timevar)){
    
    timevarnames <- sort(unique(timevar))
    correstab <- setNames(seq_along(timevarnames), timevarnames)
    
    timevarnum <- as.numeric(correstab[as.character(timevar)])
    
    K = sort(unique(timevarnum))
    
    logll3 <- function(sigma, nu, theta, time, event, covatime, hP, w, K){
      
      value = 0
      for(k in K){
        sigmak <- sigma[k]
        nuk <- nu[k]
        thetak <- theta[k]
        timek <- time[covatime == k]
        eventk <- event[covatime == k]
        hPk <- hP[covatime == k]
        wk <- w[covatime == k]
        value_strate <- -1*sum( wk*( eventk*log(hPk+ 
                                                  (1/thetak)*(1+(timek/sigmak)^nuk)^((1/thetak)-1)*(nuk/sigmak)*(timek/sigmak)^(nuk-1) )
                                     + (1-(1+(timek/sigmak)^nuk)^(1/thetak)) ) )
        value <- value +value_strate
      }
      return(value)
      
    }
    
    if (dist=="genweibull"){
      for (i in timevarnames){
        assign(paste0("param_names", i), c(paste0("log sigma_",i), paste0("log nu_",i), paste0("log theta_",i)))}
      
      init1 <- rep(0,3*length(K))
      
      loglik3 <- function(par, time, event, covatime, hP, w, K){
        
        sigma <- exp(par[1:length(K)])
        nu <- exp(par[(length(K)+1):(2*length(K))]) 
        theta <- exp(par[(2*length(K)+1):(3*length(K))])  
        return(logll3(sigma, nu, theta, time, event, covatime, hP, w, K)) }
      
      label <- covnames
      for (j in 1:3){
        for (i in timevarnames) {
          label <- c(label, get(paste0("param_names", i))[j] ) }
      }
    }
    
    if (dist=="weibull"){
      for (i in timevarnames){
        assign( paste0("param_names", i), c(paste0("log sigma_",i), paste0("log nu_",i)))}
      
      init1 <- rep(0,2*length(K))
      
      loglik3 <- function(par, time, event, covatime, hP, w, K){
        sigma <- exp(par[1:length(K)])
        nu <- exp(par[(length(K)+1):(2*length(K))])
        theta <- rep(1,length(K)) 
        return(logll3(sigma, nu, theta, time, event, covatime, hP, w,K)) }
      
      label <- covnames
      for (j in 1:2){
        for (i in timevarnames) {
          label <- c(label, get(paste0("param_names", i))[j] ) }
      }
    }
    
    if (dist=="exponential"){
      for (i in timevarnames){
        assign(paste0("param_names", i), c(paste0("log sigma_",i)))}
      
      init1 <- rep(0,length(K))
      
      loglik3 <- function(par, time, event, covatime, hP, w, K){
        
        sigma <- exp(par[1:length(K)])
        nu <- rep(1,length(K))
        theta <- rep(1,length(K))  
        return(logll3(sigma, nu, theta, time, event, covatime, hP, w, K)) }
      
      label <- covnames
      for (i in timevarnames) {
        label <- c(label, get(paste0("param_names", i)))}
      
    }
    
    suppressWarnings({
      logllmax1 <- optim(par = init1, fn = loglik3, time = time, event = event,
                         covatime = timevarnum, hP = hP, w = weights
                         , K= K)
      
      indic <- 0
      while(indic <= 5){
        ll_val <- logllmax1$value
        logllmax1 <- optim(par = logllmax1$par, fn = loglik3, time = time, 
                           event = event, covatime = timevarnum,
                           hP = hP, w = weights, K = K)
        delta <- ll_val - logllmax1$value
        if(delta ==0) {indic = indic + 1}
      }
      
      logllmax1 <- optim(par = logllmax1$par, fn = loglik3, time = time, 
                         event = event, covatime = timevarnum,
                         hP = hP, w = weights, K= K,
                         hessian = TRUE)
    })
    
  }
  
  logll0 <- function(sigma, nu, theta, time, event, hP, w){
    return(-1*sum(w*(event*log(hP+(
      (1/theta)*(1+(time/sigma)^nu)^((1/theta)-1)*(nu/sigma)*(time/sigma)^(nu-1) ) )
      + (1-(1+(time/sigma)^nu)^(1/theta)) ) ) ) }
  
  if (dist=="genweibull"){
    if(is.null(label)){
      label <- c(covnames,"log sigma", "log nu", "log theta")}
    
    init0 <- rep(0,3)
    
    loglik0 <- function(par, time, event, cova, hP, w){
      sigma <- exp(par[1]); nu <- exp(par[2]);  theta <- exp(par[3]) 
      return(logll0(sigma, nu, theta, time, event, hP, w)) }
    
    method <- "Nelder-Mead" 
    lower <- -Inf
    upper <-  Inf
  }
  
  if (dist=="weibull"){
    if(is.null(label)){
      label <- c(covnames,"log sigma", "log nu")}
    
    init0 <- rep(0,2)
    
    loglik0 <- function(par, time, event, cova, hP, w){
      sigma <- exp(par[1]); nu <- exp(par[2]);  theta <- 1
      return(logll0(sigma, nu, theta, time, event, hP, w)) }
    
    method <- "Nelder-Mead" 
    lower <- -Inf
    upper <-  Inf
  }
  
  if (dist=="exponential"){
    if(is.null(label)){
      label <- c(covnames,"log sigma")}
    
    init0 <- rep(0,1)
    
    loglik0 <- function(par, time, event, cova, hP, w){
      sigma <- exp(par[1]); nu <- 1;  theta <- 1
      return(logll0(sigma, nu, theta, time, event, hP, w)) }
    
    method <- "Brent"
    lower <- -2500
    upper <-  2500
  }
  suppressWarnings({
    
    logllmax0 <- optim(par = init0, fn = loglik0, time = time, event = event,
                       hP = hP, w = weights, method = method, lower = lower,
                       upper = upper)
    
    indic <- 0
    while(indic <= 5){
      ll_val <- logllmax0$value
      logllmax0 <- optim(par = logllmax0$par, fn = loglik0, time = time, 
                         event = event, hP = hP, w = weights, method = method, 
                         lower = lower, upper = upper)
      delta <- ll_val - logllmax0$value
      if(delta ==0) {indic = indic + 1}
    }
    
    logllmax0 <- optim(par = logllmax0$par, fn = loglik0, time = time, 
                       event = event, cova = cova, hP = hP, w = weights,
                       hessian = TRUE, method = method, lower = lower,
                       upper = upper)
  })
  
  if (!is.null(covnames)){
    
    t.table <- data.frame(coef = logllmax1$par,
                          ecoef = exp(logllmax1$par),
                          se = sqrt(diag(solve(logllmax1$hessian))),
                          z = logllmax1$par/sqrt(diag(solve(logllmax1$hessian))),
                          p = 2*(1-pnorm(abs(logllmax1$par/sqrt(diag(solve(logllmax1$hessian)))), 0, 1)),
                          row.names = label)
  }
  
  if(is.null(covnames) & is.null(timevar)){
    
    t.table <- data.frame(coef = logllmax0$par,
                          ecoef = exp(logllmax0$par),
                          se = sqrt(diag(solve(logllmax0$hessian))),
                          z = logllmax0$par/sqrt(diag(solve(logllmax0$hessian))),
                          p = 2*(1-pnorm(abs(logllmax0$par/sqrt(diag(solve(logllmax0$hessian)))), 0, 1)),
                          row.names = label)
  }
  
  if(is.null(covnames) & !is.null(timevar)){
    
    t.table <- data.frame(coef = logllmax1$par,
                          ecoef = exp(logllmax1$par),
                          se = sqrt(diag(solve(logllmax1$hessian))),
                          z = logllmax1$par/sqrt(diag(solve(logllmax1$hessian))),
                          p = 2*(1-pnorm(abs(logllmax1$par/sqrt(diag(solve(logllmax1$hessian)))), 0, 1)),
                          row.names = label)
  }
  
  names(t.table) <- c("coef", "exp(coef)", "se(coef)", "z", "p")
  
  coefficients <- t.table$coef
  names(coefficients) <- label
  
  betaestim <- coefficients[covnames]
  lp <- cova %*% betaestim
  
  dimnames(cova)[[2]] <- covnames
  
  var <- if(!is.null(covnames) || (is.null(covnames) & !is.null(timevar)) ){solve(logllmax1$hessian)
  }else{solve(logllmax0$hessian)}
  loglik <- if (!is.null(covnames) || (is.null(covnames) & !is.null(timevar)) ){c(-1*logllmax1$value, -1*logllmax0$value)
  }else{-1*logllmax0$value}
  if(length(loglik)==2){names(loglik) <- c("Model", "Null model")}else{
    names(loglik) <- c("Null Model")
  }
  
  res <- list(
    formula = formula,
    dist = dist,
    coefficients =  coefficients,
    var = var,
    t.table = t.table,
    loglik = loglik,
    linear.predictors = as.vector(lp),
    missing = !na,
    n = length(time),
    nevent = sum(event),
    y = cbind(time = time, status = event),
    x = cova,
    ays = data.frame(age = age, year = year, sex = sex),
    call = call
  )
  if (!is.null(xlevels)) {
    res$correstab <- correstab
    res$xlevels <- xlevels
    res$levelsval <- data[,strata_var]  }
  class(res) <- "survivalNET"
  return(res)
}



