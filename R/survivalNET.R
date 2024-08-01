
survivalNET <- function(formula, data, ratetable, age, year, sex, 
                        dist="weibull", strata=NULL, weights=NULL)
{
  
  call = match.call()
  
  ####### check errors
  
  if (missing(formula)) stop("a formula argument is required")
  if (missing(data)) stop("a data argument is required")
  if (missing(ratetable)) stop("a table argument is required")
  if (missing(age)) stop("an age argument is required")
  if (missing(sex)) stop("a sex argument is required")
  if (missing(year)) stop("a year argument is required")
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
  cova <- as.matrix(data[,gsub("\\+", "", attr(terms(formula), "term.labels"))])
  age <- data[,age] 
  sex <- data[,sex] 
  year <- data[,year] 
  if(is.null(weights)){weights = rep(1,dim(data[1]))}
  if(max(time)<30) stop("The event times must be expressed in days (max(", as.character(formula[[2]][2]),") is less than 30 days).")
  if(!is.numeric(age)) stop("Argument 'age' must be a numeric vector")
  if(max(age)<30) stop("The ages must be expressed in days (max(", as.character(formula[[2]][2]),") is less than 30 days).")
  if(!is.character(sex)) stop("Argument 'sex' must be a character string")
  if(min(names(table(sex)) %in% c("female","male", NA))==0) stop("Argument 'sex' must be 'male' or 'female'")
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
   
    if (length( attr( terms(formula), "term.labels")) != 0){
   
      logll1 <- function(sigma, nu, theta, beta, time, event, cova, hP, w){
     return(-1*sum( w*( event*log(hP+exp(cova%*%beta)*(
       (1/theta)*(1+(time/sigma)^nu)^((1/theta)-1)*(nu/sigma)*(time/sigma)^(nu-1) ) )
       + exp(cova%*%beta)*(1-(1+(time/sigma)^nu)^(1/theta)) ) ) )}
    }
   
    logll0 <- function(sigma, nu, theta, time, event, hP, w){
       return(-1*sum(w*(event*log(hP+(
         (1/theta)*(1+(time/sigma)^nu)^((1/theta)-1)*(nu/sigma)*(time/sigma)^(nu-1) ) )
         + (1-(1+(time/sigma)^nu)^(1/theta)) ) ) ) }
   
   if (dist=="genweibull"){
     label <- c("log sigma", "log nu", "log theta", gsub("\\+", "", attr(terms(formula), "term.labels")))
     
     if (length( attr( terms(formula), "term.labels")) != 0){
       
       init1 <- c(rep(0,3),rep(0,dim(cova)[2]))
      
       loglik1 <- function(par, time, event, cova, hP, w){
            sigma <- exp(par[1]); nu <- exp(par[2]);  theta <- exp(par[3]); beta <- par[4:(3+dim(cova)[2])] 
            return(logll1(sigma, nu, theta, beta, time, event, cova, hP, w)) }
     }
     
    
     init0 <- rep(0,3)
     
     loglik0 <- function(par, time, event, cova, hP, w){
       sigma <- exp(par[1]); nu <- exp(par[2]);  theta <- exp(par[3]) 
       return(logll0(sigma, nu, theta, time, event, hP, w)) }
    }
  
   
   if (dist=="weibull"){
     label <- c("log sigma", "log nu", gsub("\\+", "", attr(terms(formula), "term.labels")))
     
     if (length( attr( terms(formula), "term.labels")) != 0){
       
       init1 <- c(rep(0,2),rep(0,dim(cova)[2]))
       
       loglik1 <- function(par, time, event, cova, hP, w){
         sigma <- exp(par[1]); nu <- exp(par[2]);  theta <- 1; beta <- par[3:(2+dim(cova)[2])] 
         return(logll1(sigma, nu, theta, beta, time, event, cova, hP, w)) }
     }
     
     init0 <- rep(0,2)
     
     loglik0 <- function(par, time, event, cova, hP, w){
       sigma <- exp(par[1]); nu <- exp(par[2]);  theta <- 1
       return(logll0(sigma, nu, theta, time, event, hP, w)) }
   }
   
   if (dist=="exponential"){
     label <- c("log sigma", gsub("\\+", "", attr(terms(formula), "term.labels")))
     
     if (length( attr( terms(formula), "term.labels")) != 0){
       
       init1 <- c(rep(0,1),rep(0,dim(cova)[2]))
       
       loglik1 <- function(par, time, event, cova, hP, w){
         sigma <- exp(par[1]); nu <- 1;  theta <- 1; beta <- par[2:(1+dim(cova)[2])] 
         return(logll1(sigma, nu, theta, beta, time, event, cova, hP, w)) }
     
     }
     
     init0 <- rep(0,1)
     
     loglik0 <- function(par, time, event, cova, hP, w){
       sigma <- exp(par[1]); nu <- 1;  theta <- 1
       return(logll0(sigma, nu, theta, time, event, hP, w)) }
   }
   
   if (dist == "exponential") {
     method <- "Brent"
     lower <- -2500
     upper <-  2500
   } else {
     method <- "Nelder-Mead" 
     lower <- -Inf
     upper <-  Inf
   }
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
   
   if (length( attr( terms(formula), "term.labels")) != 0){
     
     
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
   }
   
   if (length( attr( terms(formula), "term.labels")) != 0){
     
     t.table <- data.frame(coef = logllmax1$par,
                           ecoef = exp(logllmax1$par),
                           se = sqrt(diag(solve(logllmax1$hessian))),
                           z = logllmax1$par/sqrt(diag(solve(logllmax1$hessian))),
                           p = 2*(1-pnorm(abs(logllmax1$par/sqrt(diag(solve(logllmax1$hessian)))), 0, 1)),
                           row.names = label)
   }
   else{
     
     t.table <- data.frame(coef = logllmax0$par,
                           ecoef = exp(logllmax0$par),
                           se = sqrt(diag(solve(logllmax0$hessian))),
                           z = logllmax0$par/sqrt(diag(solve(logllmax0$hessian))),
                           p = 2*(1-pnorm(abs(logllmax0$par/sqrt(diag(solve(logllmax0$hessian)))), 0, 1)),
                           row.names = label)
   }
  
  names(t.table) <- c("coef", "exp(coef)", "se(coef)", "z", "p")
  
  coefficients <- t.table$coef
  names(coefficients) <- label
  
  betaestim <- coefficients[names(coefficients) %in% gsub("\\+", "", attr(terms(formula), "term.labels"))]
  lp <- cova %*% betaestim
  
  dimnames(cova)[[2]] <- gsub("\\+", "", attr(terms(formula), "term.labels"))
  
  var = if (length( attr( terms(formula), "term.labels")) != 0){solve(logllmax1$hessian)
  }else{solve(logllmax0$hessian)}
  loglik = if (length( attr( terms(formula), "term.labels")) != 0){c(-1*logllmax1$value, -1*logllmax0$value)
  }else{-1*logllmax0$value}
    
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
    asy = data.frame(age = age, sex = sex, year = year),
    call = call
  )
  class(res) <- "survivalNET"
  return(res)
  }



