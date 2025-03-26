
survivalFLEXNET <- function(formula, data, ratetable, m=3, mpos = NULL, mquant = NULL, init = NULL, 
                            delta_th = 0, weights=NULL)
{
  
  ####### check errors
  
  if (missing(formula)) stop("a formula argument is required")
  if (missing(data)) stop("a data argument is required")
  if (missing(ratetable)) stop("a ratetable argument is required")
  if (as.character(class(formula)) != "formula") stop("The first argument must be a formula")
  if (as.character(class(data)) != "data.frame") stop("The second argument must be a data frame")
  if (length(dim(ratetable))!=3) stop("The life table must have 3 dimensions: age, year, sex")
  if (dim(ratetable)[3]!=2) stop("The life table must have 3 dimensions: age, year, sex")
  if(!is.null(init)){if(!is.numeric(init))stop("Argument 'init' must be a vector of numeric values.") }
  if(!is.numeric(delta_th))stop("'delta_th' must be numeric.")
  if(length(delta_th) != 1) stop("'delta_th' must be a single value.") 
  if(!is.null(mpos) & !is.null(mquant))warning("'mpos' and 'mquant' have both been specified. 'mpos' values have been chosen over 'mquant' quantiles.") #ligne : mquant = NULL
  
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
    return(list(age = age, year = year, sex = sex))
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
  if(min(names(table(sex)) %in% c("female","male", NA))==0) stop("'sex' must be 'male' or 'female'")
  # if(!is.date(year)) stop("The values for 'year' must be of the 'date' class")
  if(!(is.numeric(year) || inherits(year, "Date") || inherits(year, "POSIXct") || inherits(year, "POSIXlt")))stop("'year' must be of class numeric, Date, POSIXct, or POSIXlt.")
  if(!is.null(weights)){
    if(!is.numeric(weights)) stop("Argument 'weights' must be a numeric vector")
    if(length(weights)!=dim(data)[1]) stop("Argument 'weights' must have the same length as the number of rows of the 'data' argument. (", dim(data)[1],")")}
  if(!is.null(mpos) & !is.null(mquant)){mquant = NULL}
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
    
    logll1 <- function(beta, gamma, time, event, cova, hP, w, m, mpos, mquant){
      return(-1*sum( w*(
        event * log(hP + (1/time)*splinecubeP(time, gamma, m, mpos, mquant)$spln *
                      exp(splinecube(time, gamma, m, mpos, mquant)$spln + cova %*% beta) ) -
          exp(splinecube(time, gamma, m, mpos, mquant)$spln + cova %*% beta)
      ) ) ) } ###p5489 StatinMed P.Nelson PC. Lambert Flexible Parametric models for relative survival
    
    
    if(m == 0){gamma_names = NULL
    }else{gamma_names <- paste0("gamma", 2:(m+1))}
    
    label <- c(covnames, "gamma0",
                   "gamma1", gamma_names)
    
    if(!is.null(init)){
      if(length(init) != dim(cova)[2]+m+2) stop("'init' length must be ", dim(cova)[2]+m+2,
                                      " ( ", dim(cova)[2]," covariate(s) and ",m+2," parameters for the Restricted Cubic Spline).")
      init1 <- init
    }else{init1 <- c(rep(0,dim(cova)[2]+m+2))}
    
    loglik1 <- function(par, time, event, cova, hP, w, m, mpos, mquant){
      beta <- par[1:dim(cova)[2]]
      gamma <- par[(dim(cova)[2]+1):length(par)]
      return(logll1(beta, gamma, time, event, cova, hP, w, m, mpos, mquant)) }
    
    suppressWarnings({
    logllmax1 <- optim(par = init1, fn = loglik1, time = time, event = event,
                       cova = cova, hP = hP, w = weights, m = m, mpos = mpos, mquant = mquant)
    
    indic <- 0
    while(indic <= 5){
      ll_val <- logllmax1$value
      logllmax1 <- optim(par = logllmax1$par, fn = loglik1, time = time, 
                         event = event, cova = cova, hP = hP, w = weights,
                         m = m, mpos = mpos, mquant = mquant)
      delta <- ll_val - logllmax1$value
      if(delta_th == 0){
        if(delta == delta_th) {indic = indic + 1}
      }else{ 
        if(0 < delta & delta <= delta_th) {indic = indic + 1}
      }      }
    
    logllmax1 <- optim(par = logllmax1$par, fn = loglik1, time = time, 
                       event = event, cova = cova, hP = hP, w = weights,
                       m = m, mpos = mpos, mquant = mquant,
                       hessian = TRUE)
    })
  
  }
  
  if(!is.null(covnames) & !is.null(timevar)){
    
    timevarnames <- sort(unique(timevar))
    correstab <- setNames(seq_along(timevarnames), timevarnames)
    
    timevarnum <- as.numeric(correstab[as.character(timevar)])
    
    K = sort(unique(timevarnum))
    
    logll2 <- function(beta, gamma, time, event, cova, covatime, hP, w, m, mpos, mquant, K){
      
      value = 0
      for(k in K){
        gammak <- gamma[,k] 
        timek <- time[covatime == k]
        eventk <- event[covatime == k]
        hPk <- hP[covatime == k]
        covak <- cova[covatime == k]
        wk <- w[covatime == k]
        value_strate <- -1*sum(wk * (eventk * log(hPk + (1/timek)*splinecubeP(timek, gammak, m, mpos, mquant)$spln*
                                     exp(splinecube(timek, gammak, m, mpos, mquant)$spln + as.matrix(covak) %*% beta) ) -
          exp(splinecube(timek, gammak, m, mpos, mquant)$spln + as.matrix(covak) %*% beta)
        )
        )
        value <- value +value_strate
      }
      return(value)
    }
    
    if (m == 0) {
      for (i in timevarnames) {
        assign(paste0("gamma_names", i), NULL)
      }
    }else {
      for (i in timevarnames) {
        assign(paste0("gamma_names", i), paste0("gamma", i, "_", 2:(m+1)))
      }
    }
    
    label <- covnames
    for (i in timevarnames) {
      label <- c(label, paste0("gamma", i, "_0"), paste0("gamma", i, "_1"), get(paste0("gamma_names", i)))
    }
    
    if(!is.null(init)){
      if(length(init) != dim(cova)[2]+(length(K)*(m+2))) stop("'init' length must be ", dim(cova)[2]+m+2,
                                                " ( ", dim(cova)[2]," covariate(s) and ",(length(K)*(m+2))," parameters for the Restricted Cubic Spline for each level (",length(K),") of the covariate with a time dependant effect).")
      init1 <- init
    }else{init1 <- c(rep(0,dim(cova)[2]+(length(K)*(m+2)))) }
    
    loglik2 <- function(par, time, event, cova, covatime, hP, w, m, mpos, mquant, K){
      beta <- par[1:dim(cova)[2]]
      gamma <-matrix(par[(dim(cova)[2]+1):(dim(cova)[2]+(length(K)*(m+2)))], ncol = length(K))
      
      return(logll2(beta, gamma, time, event, cova, covatime, hP, w, m, mpos, mquant, K)) }
    
    suppressWarnings({
    logllmax1 <- optim(par = init1, fn = loglik2, time = time, event = event,
                       cova = cova, covatime = timevarnum, hP = hP, w = weights
                       , m = m, mpos = mpos, mquant = mquant, K= K)
    
    indic <- 0
    while(indic <= 5){
      last_logllmax1 <- logllmax1
      ll_val <- logllmax1$value
      logllmax1 <- optim(par = logllmax1$par, fn = loglik2, time = time, 
                         event = event, cova = cova, covatime = timevarnum,
                         hP = hP, w = weights, m = m, mpos = mpos, mquant = mquant, K = K)
      delta <- ll_val - logllmax1$value
      if(delta_th == 0){
        if(delta == delta_th) {indic = indic + 1}
      }else{ 
        if(0 < delta & delta <= delta_th) {indic = indic + 1}
      }      }
    
    tryCatch({logllmax1 <- optim(par = logllmax1$par, fn = loglik2, time = time, 
                       event = event, cova = cova, covatime = timevarnum,
                       hP = hP, w = weights, m = m, mpos = mpos, mquant = mquant, K= K,
                       hessian = TRUE)
             
    }, error = function(e){logllmax1 <- last_logllmax1})
    })
  }
  
  if (is.null(covnames) & !is.null(timevar)){
    
    timevarnames <- sort(unique(timevar))
    correstab <- setNames(seq_along(timevarnames), timevarnames)
    
    timevarnum <- as.numeric(correstab[as.character(timevar)])
    
    K = sort(unique(timevarnum))
    
    logll3 <- function(gamma, time, event, covatime, hP, w, m, mpos, mquant, K){
      
      value = 0
      for(k in K){
        gammak <- gamma[,k] 
        timek <- time[covatime == k]
        eventk <- event[covatime == k]
        hPk <- hP[covatime == k]
        wk <- w[covatime == k]
        value_strate <- -1*sum(wk * (eventk * log(hPk + (1/timek)*splinecubeP(timek, gammak, m, mpos, mquant)$spln*
                                                    exp(splinecube(timek, gammak, m, mpos, mquant)$spln ) ) -
                                       exp(splinecube(timek, gammak, m, mpos, mquant)$spln)
        )
        )
        value <- value +value_strate
      }
      return(value)
    }
    
    if (m == 0) {
      for (i in timevarnames) {
        assign(paste0("gamma_names", i), NULL)
      }
    }else {
      for (i in timevarnames) {
        assign(paste0("gamma_names", i), paste0("gamma", i, "_", 2:(m+1)))
      }
    }
    
    label <- covnames
    for (i in timevarnames) {
      label <- c(label, paste0("gamma", i, "_0"), paste0("gamma", i, "_1"), get(paste0("gamma_names", i)))
    }
    
    if(!is.null(init)){
      if(length(init) != (length(K)*(m+2))) stop("'init' length must be ",(length(K)*(m+2))," (",(length(K)*(m+2))," parameters for the Restricted Cubic Spline for each level (",length(K),") of the covariate with a time dependant effect).")
      init1 <- init
    }else{init1 <- rep(0,(length(K)*(m+2)))}
    
    loglik3 <- function(par, time, event, covatime, hP, w, m, mpos, mquant, K){
     
       gamma <-matrix(par[1:((length(K)*(m+2)))], ncol = length(K))
      
      return(logll3(gamma, time, event, covatime, hP, w, m, mpos, mquant, K)) }
    
    suppressWarnings({
      logllmax1 <- optim(par = init1, fn = loglik3, time = time, event = event,
                          covatime = timevarnum, hP = hP, w = weights
                         , m = m, mpos = mpos, mquant = mquant, K= K)
      
      indic <- 0
      while(indic <= 5){
        last_logllmax1 <- logllmax1
        ll_val <- logllmax1$value
        logllmax1 <- optim(par = logllmax1$par, fn = loglik3, time = time, 
                           event = event, covatime = timevarnum,
                           hP = hP, w = weights, m = m, mpos = mpos, mquant = mquant, K = K)
        delta <- ll_val - logllmax1$value
        if(delta_th == 0){
          if(delta == delta_th) {indic = indic + 1}
        }else{ 
          if(0 < delta & delta <= delta_th) {indic = indic + 1}
        }        }
      
      tryCatch({logllmax1 <- optim(par = logllmax1$par, fn = loglik3, time = time, 
                         event = event, covatime = timevarnum,
                         hP = hP, w = weights, m = m, mpos = mpos, mquant = mquant, K= K,
                         hessian = TRUE)}
      , error = function(e){logllmax1 <- last_logllmax1})
    })
  }

  
  #NULL model
  logll0 <- function(gamma, time, event, hP, w, m, mpos, mquant){
    return(-1*sum(w*(
      event * log(hP + (1/time)*splinecubeP(time, gamma, m, mpos, mquant)$spln *
                    exp(splinecube(time, gamma, m, mpos, mquant)$spln) ) -
        exp(splinecube(time, gamma, m, mpos, mquant)$spln)
     ) )
    ) 
    }
  
  if(m == 0){gamma_names = NULL
  }else{gamma_names <- paste0("gamma", 2:(m+1))}
  
  labelNULL <- c(covnames, "gamma0",
             "gamma1", gamma_names)
  
  if(is.null(covnames) & is.null(timevar)){
    if(!is.null(init)){
      if(length(init) != (m+2)) stop("'init' length must be ",(m+2)," (",(m+2)," parameters for the Restricted Cubic Spline).")
    init0 <- init}else{init0 <- rep(0,m+2)
    }}else{
      init0 <- rep(0,m+2)
    }
  
  
  loglik0 <- function(par, time, event, hP, w, m, mpos, mquant){
    gamma <- par
    return(logll0(gamma, time, event, hP, w, m, mpos, mquant)) }
  
  suppressWarnings({
    
  logllmax0 <- optim(par = init0, fn = loglik0, time = time, event = event,
                     hP = hP, w = weights, m = m, mpos = mpos, mquant = mquant)
  
  indic <- 0
  while(indic <= 5){
    ll_val <- logllmax0$value
    logllmax0 <- optim(par = logllmax0$par, fn = loglik0, time = time, 
                       event = event, hP = hP, w = weights, m = m, mpos = mpos, mquant = mquant)
    delta <- ll_val - logllmax0$value
    if(delta_th == 0){
      if(delta == delta_th) {indic = indic + 1}
    }else{ 
      if(0 < delta & delta <= delta_th) {indic = indic + 1}
    }  
    }
  
  logllmax0 <- optim(par = logllmax0$par, fn = loglik0, time = time, 
                     event = event, hP = hP, w = weights, 
                     hessian = TRUE, m = m, mpos = mpos, mquant = mquant)
  })
  
  ##récupération de mpos et mquant si ils n'ont pas été spécifiés dans la formule.
  
  if(is.null(mpos)){
    if(is.null(mquant)){
      a <- c()
      for(i in (0:(m+1))){
        a <- c(a,i/(m+1))}
      mpos <- quantile(log(time), probs = a)
      mpos <- as.numeric(mpos)
      mquant <- a 
    }else{
      a <- c(mquant)
      mpos <- quantile(log(time), probs = a)
    }
  }
  
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
                          row.names = labelNULL)
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
  loglik <- if(!is.null(covnames) || (is.null(covnames) & !is.null(timevar)) ){c(-1*logllmax1$value, -1*logllmax0$value)
  }else{-1*logllmax0$value}
  if(length(loglik)==2){names(loglik) <- c("Model", "Null model")}else{
    names(loglik) <- c("Null Model")
  }
  
  res <- list(
    formula = formula,
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
    m = m,
    mpos = mpos,
    mquant = mquant
  )
  if (!is.null(xlevels)) {
    res$correstab <- correstab
    res$xlevels <- xlevels
    res$levelsval <- data[,strata_var]
  }
  class(res) <- "survivalNET"
  return(res)
}






