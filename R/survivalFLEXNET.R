
survivalFLEXNET <- function(formula, data, ratetable, age, year, sex,
                        m=3, mpos = NULL, strata=NULL, weights=NULL)
{
  
  ####### check errors
  
  if (missing(formula)) stop("a formula argument is required")
  if (missing(data)) stop("a data argument is required")
  if (missing(ratetable)) stop("a table argument is required")
  if (missing(age)) stop("an age argument is required")
  if (missing(sex)) stop("a sex argument is required")
  if (missing(year)) stop("a year argument is required")
  if (as.character(class(formula)) != "formula") stop("The first argument must be a formula")
  if (as.character(class(data)) != "data.frame") stop("The second argument must be a data frame")
  if (length(dim(ratetable))!=3) stop("The life table must have 3 dimensions: age, year, sex")
  if (dim(ratetable)[3]!=2) stop("The life table must have 3 dimensions: age, year, sex")

  ####### data management
  
  time <- data[,as.character(formula[[2]][2])] # Thomas: test initial une unite en jours comme relsurv
  event <- data[,as.character(formula[[2]][3])]
  cova <- as.matrix(data[,gsub("\\+", "", attr(terms(formula), "term.labels"))])
  age <- data[,age] # Thomas: test initial une unite en jours 
  sex <- data[,sex] 
  year <- data[,year]
  if(max(time)<30) stop("The event times must be expressed in days (max(", as.character(formula[[2]][2]),") is less than 30 days).")
  if(!is.numeric(age)) stop("Argument 'age' must be a numeric vector")
  if(max(age)<30) stop("The ages must be expressed in days (max(", as.character(formula[[2]][2]),") is less than 30 days).")
  if(!is.character(sex)) stop("Argument 'sex' must be a character string")
  if(min(names(table(sex)) %in% c("female","male", NA))==0) stop("Argument 'sex' must be 'male' or 'female'")
  if(!is.date(year)) stop("The values for 'year' must be of the 'date' class")
  for(i in colnames(cova)){ if(!is.numeric(cova[,i])) stop("All covariates must be numeric")  }
  
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
   
   
  logll1 <- function(beta, gamma, time, event, cova, hP, m, mpos){
   return(-1*sum(
     event * log(hP + (1/time)*splinecubeP(time, gamma, m, mpos)$spln *
                 exp(splinecube(time, gamma, m)$spln + cova %*% beta) ) -
       exp(splinecube(time, gamma, m)$spln + cova %*% beta)
   )) }
   
  logll0 <- function(gamma, time, event, cova, hP, m, mpos){
     return(-1*sum(
       event * log(hP + (1/time)*splinecubeP(time, gamma, m, mpos)$spln *
                     exp(splinecube(time, gamma, m)$spln) ) -
         exp(splinecube(time, gamma, m)$spln)
     )) }
   
  if(m == 0){gamma_names = NULL
  }else{gamma_names <- paste0("gamma", 2:(m+1))}
   
  label <- c(gsub("\\+", "", attr(terms(formula), "term.labels")), "gamma0",
               "gamma1", gamma_names)
   
  init1 <- c(rep(0,dim(cova)[2]+m+2))
     
  loglik1 <- function(par, time, event, cova, hP, m, mpos){
       beta <- par[1:dim(cova)[2]]
       gamma <- par[(dim(cova)[2]+1):length(par)]
       return(logll1(beta, gamma, time, event, cova, hP, m, mpos)) }
     
  init0 <- rep(0,m+2)
     
  loglik0 <- function(par, time, event, cova, hP, m, mpos){
       gamma <- par
       return(logll0(gamma, time, event, cova, hP, m, mpos)) }

   
  logllmax0 <- optim(par = init0, fn = loglik0, time = time, event = event,
                      hP = hP, m = m, mpos = mpos)
   
  indic <- 0
  while(indic <= 5){
     ll_val <- logllmax0$value
     logllmax0 <- optim(par = logllmax0$par, fn = loglik0, time = time, 
                        event = event, hP = hP, m = m, mpos = mpos)
     delta <- ll_val - logllmax0$value
     if(delta ==0) {indic = indic + 1}
  }
   
  logllmax0 <- optim(par = logllmax0$par, fn = loglik0, time = time, 
                      event = event, cova = cova, hP = hP, hessian = TRUE,
                      m = m, mpos = mpos)
   
  logllmax1 <- optim(par = init1, fn = loglik1, time = time, event = event,
                     cova = cova, hP = hP, m = m, mpos = mpos)
   
  indic <- 0
  while(indic <= 5){
    ll_val <- logllmax1$value
    logllmax1 <- optim(par = logllmax1$par, fn = loglik1, time = time, 
                       event = event, cova = cova, hP = hP, m = m, mpos = mpos)
    delta <- ll_val - logllmax1$value
    if(delta ==0) {indic = indic + 1}
  }
  
  logllmax1 <- optim(par = logllmax1$par, fn = loglik1, time = time, 
                     event = event, cova = cova, hP = hP, m = m, mpos = mpos,
                     hessian = TRUE)
  
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
  
  dimnames(cova)[[2]] <- gsub("\\+", "", attr(terms(formula), "term.labels"))
  
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
  class(res) <- "survivalFLEXNET"
  return(res)
  }



