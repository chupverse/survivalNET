
survivalFLEXNET <- function(formula, data, ratetable, m=3, mpos = NULL, 
                            weights=NULL)
{
  
  ####### check errors
  
  if (missing(formula)) stop("a formula argument is required")
  if (missing(data)) stop("a data argument is required")
  if (missing(ratetable)) stop("a table argument is required")
  if (as.character(class(formula)) != "formula") stop("The first argument must be a formula")
  if (as.character(class(data)) != "data.frame") stop("The second argument must be a data frame")
  if (length(dim(ratetable))!=3) stop("The life table must have 3 dimensions: age, year, sex")
  if (dim(ratetable)[3]!=2) stop("The life table must have 3 dimensions: age, year, sex")
  
  ####### data management
  
  time <- data[,as.character(formula[[2]][2])] # Thomas: test initial une unite en jours comme relsurv
  event <- data[,as.character(formula[[2]][3])]
  all_terms <- attr(terms(formula), "term.labels")
  strata_terms <- grep("strata\\(", all_terms, value = TRUE)
  ratetable_terms <- grep("ratetable\\(", all_terms, value = TRUE)
  if(length(ratetable_terms) == 0) stop("Error: The formula must contain a ratetable() term.")
  covnames <- setdiff(all_terms, c(strata_terms, ratetable_terms))
  cova <- as.matrix(data[,covnames])
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
  if(is.null(unlist(lapply(strata_terms, extract_vars)))){
    timevar = unlist(lapply(strata_terms, extract_vars))
  }
  if(!is.null(unlist(lapply(strata_terms, extract_vars)))){
  timevar <- data[,unlist(lapply(strata_terms, extract_vars))]}
  if(!is.null(timevar)){
    if(length(unique(timevar))>15) stop("The variable with a time-dependant effect has too many categories (>15)")
  }
  ratetable_vars <- assign_ratetable_vars(unlist(lapply(ratetable_terms, extract_vars_within_function)))
  age <- data[,ratetable_vars$age] # Thomas: test initial une unite en jours 
  year <- data[,ratetable_vars$year]
  sex <- data[,ratetable_vars$sex] 
  if(is.null(weights)){weights = rep(1,dim(data)[1])}
  if(max(time)<30) stop("The event times must be expressed in days (max(", as.character(formula[[2]][2]),") is less than 30 days).")
  if(!is.numeric(age)) stop("'age' must be numeric")
  if(max(age)<30) stop("The ages must be expressed in days (max(", as.character(formula[[2]][2]),") is less than 30 days).")
  if(!is.character(sex)) stop("'sex' must be a character string")
  if(min(names(table(sex)) %in% c("female","male", NA))==0) stop("'sex' must be 'male' or 'female'")
  if(!is.date(year)) stop("The values for 'year' must be of the 'date' class")
  for(i in colnames(cova)){ if(!is.numeric(cova[,i])) stop("All covariates must be numeric")  }
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
  
  #hP <- expectedhaz(ratetable, age=age, sex=sex, year=year, time=time) # compute instantaneous hazards
  
  hP <- sapply(seq_along(time), function(i) {
    expectedhaz(ratetable = ratetable, age = age[i], sex = sex[i], year = year[i], time = time[i])
  })
  
  ###### log likelihood functions
  
  if(is.null(timevar)){
    
    logll1 <- function(beta, gamma, time, event, cova, hP, w, m, mpos){
      return(-1*sum( w*(
        event * log(hP + (1/time)*splinecubeP(time, gamma, m, mpos)$spln *
                      exp(splinecube(time, gamma, m)$spln + cova %*% beta) ) -
          exp(splinecube(time, gamma, m)$spln + cova %*% beta)
      ) ) ) }
    
    if(m == 0){gamma_names = NULL
    }else{gamma_names <- paste0("gamma", 2:(m+1))}
    
    label <- c(covnames, "gamma0",
               "gamma1", gamma_names)
    
    init1 <- c(rep(0,dim(cova)[2]+m+2))
    
    loglik1 <- function(par, time, event, cova, hP, w, m, mpos){
      beta <- par[1:dim(cova)[2]]
      gamma <- par[(dim(cova)[2]+1):length(par)]
      return(logll1(beta, gamma, time, event, cova, hP, w, m, mpos)) }
    
    
    logllmax1 <- optim(par = init1, fn = loglik1, time = time, event = event,
                       cova = cova, hP = hP, w = weights, m = m, mpos = mpos)
    
    indic <- 0
    while(indic <= 5){
      ll_val <- logllmax1$value
      logllmax1 <- optim(par = logllmax1$par, fn = loglik1, time = time, 
                         event = event, cova = cova, hP = hP, w = weights,
                         m = m, mpos = mpos)
      delta <- ll_val - logllmax1$value
      if(delta ==0) {indic = indic + 1}
    }
    
    logllmax1 <- optim(par = logllmax1$par, fn = loglik1, time = time, 
                       event = event, cova = cova, hP = hP, w = weights,
                       m = m, mpos = mpos,
                       hessian = TRUE)
  }
  
  if(!is.null(timevar)){
    
    timevarnames <- sort(unique(timevar))
    correstab <- setNames(seq_along(timevarnames), timevarnames)
    
    timevarnum <- as.numeric(correstab[as.character(timevar)])
    
    K = sort(unique(timevarnum))
    
    logll1 <- function(beta, gamma, time, event, cova, covatime, hP, w, m, mpos, K){
      
      value = 0
      for(k in K){
        gammak <- gamma[,k] 
        timek <- time[covatime == k]
        eventk <- event[covatime == k]
        hPk <- hP[covatime == k]
        covak <- cova[covatime == k]
        wk <- w[covatime == k]
        value_strate <- -1*sum(wk * (eventk * log(hPk + (1/timek)*splinecubeP(timek, gammak, m, mpos)$spln*
                                     exp(splinecube(timek, gammak, m)$spln + as.matrix(covak) %*% beta) ) -
          exp(splinecube(timek, gammak, m)$spln + as.matrix(covak) %*% beta)
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
    
    init1 <- c(rep(0,dim(cova)[2]+(length(K)*(m+2))))
    
    loglik1 <- function(par, time, event, cova, covatime, hP, w, m, mpos, K){
      beta <- par[1:dim(cova)[2]]
      gamma <-matrix(par[(dim(cova)[2]+1):(dim(cova)[2]+(length(K)*(m+2)))], ncol = length(K))
      
      return(logll1(beta, gamma, time, event, cova, covatime, hP, w, m, mpos, K)) }

    logllmax1 <- optim(par = init1, fn = loglik1, time = time, event = event,
                       cova = cova, covatime = timevarnum, hP = hP, w = weights
                       , m = m, mpos = mpos, K= K)
    
    indic <- 0
    while(indic <= 5){
      ll_val <- logllmax1$value
      logllmax1 <- optim(par = logllmax1$par, fn = loglik1, time = time, 
                         event = event, cova = cova, covatime = timevarnum,
                         hP = hP, w = weights, m = m, mpos = mpos, K = K)
      delta <- ll_val - logllmax1$value
      if(delta ==0) {indic = indic + 1}
    }
    
    logllmax1 <- optim(par = logllmax1$par, fn = loglik1, time = time, 
                       event = event, cova = cova, covatime = timevarnum,
                       hP = hP, w = weights, m = m, mpos = mpos, K= K,
                       hessian = TRUE)
    
  }
  
  #NULL model
  logll0 <- function(gamma, time, event, cova, hP, w, m, mpos){
    return(-1*sum(w*(
      event * log(hP + (1/time)*splinecubeP(time, gamma, m, mpos)$spln *
                    exp(splinecube(time, gamma, m)$spln) ) -
        exp(splinecube(time, gamma, m)$spln)
     ) )
    ) 
    }
  
  init0 <- rep(0,m+2)
  
  loglik0 <- function(par, time, event, cova, hP, w, m, mpos){
    gamma <- par
    return(logll0(gamma, time, event, cova, hP, w, m, mpos)) }
  
  
  logllmax0 <- optim(par = init0, fn = loglik0, time = time, event = event,
                     hP = hP, w = weights, m = m, mpos = mpos)
  
  indic <- 0
  while(indic <= 5){
    ll_val <- logllmax0$value
    logllmax0 <- optim(par = logllmax0$par, fn = loglik0, time = time, 
                       event = event, hP = hP, w = weights, m = m, mpos = mpos)
    delta <- ll_val - logllmax0$value
    if(delta ==0) {indic = indic + 1}
  }
  
  logllmax0 <- optim(par = logllmax0$par, fn = loglik0, time = time, 
                     event = event, cova = cova, hP = hP, w = weights, 
                     hessian = TRUE, m = m, mpos = mpos)
  
  
  t.table <- data.frame(coef = logllmax1$par,
                        ecoef = exp(logllmax1$par),
                        se = sqrt(diag(solve(logllmax1$hessian))),
                        z = logllmax1$par/sqrt(diag(solve(logllmax1$hessian))),
                        p = 2*(1-pnorm(abs(logllmax1$par/sqrt(diag(solve(logllmax1$hessian)))), 0, 1)),
                        row.names = label)
  
  names(t.table) <- c("coef", "exp(coef)", "se(coef)", "z", "p")
  
  coefficients <- t.table$coef
  names(coefficients) <- label
  
  betaestim <- coefficients[names(coefficients) %in% gsub("\\+", "", attr(terms(formula), "term.labels"))]
  lp <- cova %*% betaestim
  
  dimnames(cova)[[2]] <- covnames
  
  res <- list(
    formula = formula,
    dist = dist,
    coefficients =  coefficients,
    var = solve(logllmax1$hessian),
    t.table = t.table,
    loglik = c(-1*logllmax1$value, -1*logllmax0$value),
    linear.predictors = as.vector(lp),
    missing = !na,
    n = length(time),
    nevent = sum(event),
    y = cbind(time = time, status = event),
    x = cova,
    asy = data.frame(age = age, sex = sex, year = year),
    m = m,
    mpos = mpos
  )
  class(res) <- "survivalNET"
  return(res)
}





