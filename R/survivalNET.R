
survivalNET <- function(formula, data, ratetable, age, year, sex, dist="weibull", strata=NULL, weights=NULL)
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
  if(!(dist %in% c("exponential","weibull","genweibull")))  stop("Argument 'dist' must be 'exponential', 'weibull', or 'genweibull' ")
  
  ####### data management
  
  time <- data[,as.character(formula[[2]][2])] # Thomas: test initial une unite en jours comme relsurv
  event <- data[,as.character(formula[[2]][3])]
  cova <- as.matrix(data[,gsub("\\+", "", attr(terms(formula), "term.labels"))])
  age <- data[,age] # Thomas: test initial une unite en jours 
  sex <- data[,sex] # Thomas: test initial si caractere "male" ou "female" -> rien d'autre (sauf eventuellement NA car supprimer ensuite)
  year <- data[,year] # Thomas: test initial si nombre de jours depuis 1960  
  if(!is.numeric(age)) stop("Argument 'age' must be a numeric vector")
  if(!is.character(sex)) stop("Argument 'sex' must be a character string")
  if(min(names(table(sex)) %in% c("female","male"))==0) stop("Argument 'p.sex' must be 'male' or 'female'")
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
   
   hP <- rep(-99, length(time))
   for (i in 1:length(time))
   {
     hP[i] <- expectedhaz(ratetable=ratetable, age=age[i], sex=sex[i], year=year[i], time=time[i])
   }
   
   # Thomas -> apply pour rapidite
   
   ###### log likelihood functions
   
   logll1 <- function(sigma, nu, theta, beta, time, event, cova, hP){
   return(-1*sum(event*log(hP+exp(cova%*%beta)*(
     (1/theta)*(1+(time/sigma)^nu)^((1/theta)-1)*(nu/sigma)*(time/sigma)^(nu-1) ) )
     + exp(cova%*%beta)*(1-(1+(time/sigma)^nu)^(1/theta)))) }
   
   logll0 <- function(sigma, nu, theta, time, event, hP){
     return(-1*sum(event*log(hP+(
       (1/theta)*(1+(time/sigma)^nu)^((1/theta)-1)*(nu/sigma)*(time/sigma)^(nu-1) ) )
       + (1-(1+(time/sigma)^nu)^(1/theta)))) }
   
   if (dist=="genweibull"){
     label <- c("log sigma", "log nu", "log theta", gsub("\\+", "", attr(terms(formula), "term.labels")))
     
     init1 <- c(rep(0,3),rep(0,dim(cova)[2]))
    
     loglik1 <- function(par, time, event, cova, hP){
          sigma <- exp(par[1]); nu <- exp(par[2]);  theta <- exp(par[3]); beta <- par[4:(3+dim(cova)[2])] 
          return(logll1(sigma, nu, theta, beta, time, event, cova, hP)) }
     
     init0 <- rep(0,3)
     
     loglik0 <- function(par, time, event, cova, hP){
       sigma <- exp(par[1]); nu <- exp(par[2]);  theta <- exp(par[3]) 
       return(logll0(sigma, nu, theta, time, event, hP)) }
    }
  
   
   if (dist=="weibull"){
     label <- c("log sigma", "log nu", gsub("\\+", "", attr(terms(formula), "term.labels")))
     
     init1 <- c(rep(0,2),rep(0,dim(cova)[2]))
     
     loglik1 <- function(par, time, event, cova, hP){
       sigma <- exp(par[1]); nu <- exp(par[2]);  theta <- 1; beta <- par[3:(2+dim(cova)[2])] 
       return(logll1(sigma, nu, theta, beta, time, event, cova, hP)) }
     
     init0 <- rep(0,2)
     
     loglik0 <- function(par, time, event, cova, hP){
       sigma <- exp(par[1]); nu <- exp(par[2]);  theta <- 1
       return(logll0(sigma, nu, theta, time, event, hP)) }
   }
   
   if (dist=="exponential"){
     label <- c("log sigma", gsub("\\+", "", attr(terms(formula), "term.labels")))
     
     init1 <- c(rep(0,1),rep(0,dim(cova)[2]))
     
     loglik1 <- function(par, time, event, cova, hP){
       sigma <- exp(par[1]); nu <- 1;  theta <- 1; beta <- par[2:(1+dim(cova)[2])] 
       return(logll1(sigma, nu, theta, beta, time, event, cova, hP)) }
     
     init0 <- rep(0,1)
     
     loglik0 <- function(par, time, event, cova, hP){
       sigma <- exp(par[1]); nu <- 1;  theta <- 1
       return(logll0(sigma, nu, theta, time, event, hP)) }
   }
   
   if (dist == "exponential") {
     method <- "Brent"
     lower <- -100000
     upper <-  100000
   } else {
     method <- "Nelder-Mead" 
     lower <- -Inf
     upper <-  Inf
   }
   logllmax0 <- optim(par = init0, fn = loglik0, time = time, event = event,
                      hP = hP, method = method, lower = lower, upper = upper)
   
   indic <- 0
   while(indic <= 5){
     ll_val <- logllmax0$value
     logllmax0 <- optim(par = logllmax0$par, fn = loglik0, time = time, 
                        event = event, hP = hP, method = method, 
                        lower = lower, upper = upper)
     delta <- ll_val - logllmax0$value
     if(delta ==0) {indic = indic + 1}
   }
   
   logllmax0 <- optim(par = logllmax0$par, fn = loglik0, time = time, 
                      event = event, cova = cova, hP = hP, hessian = TRUE,
                      method = method, lower = lower, upper = upper)
   
   logllmax1 <- optim(par = init1, fn = loglik1, time = time, event = event,
                     cova = cova, hP = hP)
   
  indic <- 0
  while(indic <= 5){
    ll_val <- logllmax1$value
    logllmax1 <- optim(par = logllmax1$par, fn = loglik1, time = time, event = event,
                cova = cova, hP = hP)
    delta <- ll_val - logllmax1$value
    if(delta ==0) {indic = indic + 1}
  }
  
  logllmax1 <- optim(par = logllmax1$par, fn = loglik1, time = time, event = event,
                    cova = cova, hP = hP, hessian = TRUE)
  
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
    asy = data.frame(age = age, sex = sex, year = year)
  )
  class(res) <- "survivalNET"
  return(res)
  }



