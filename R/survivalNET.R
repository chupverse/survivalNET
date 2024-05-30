
survivalNET <- function(formula, data, ratetable, age, sex, year, dist="weibull", strata=NULL, weights=NULL)
{
  
  ####### check errors
  
  if (missing(formula)) stop("a formula argument is required")
  if (missing(data)) stop("a data argument is required")
  if (missing(ratetable)) stop("a table argument is required")
  if (missing(age)) stop("an age argument is required")
  if (missing(sex)) stop("a sex argument is required")
  if (missing(year)) stop("a year argument is required")
  if (class(formula) != "formula") stop("The first argument must be a formula")
  if (class(data) != "data.frame") stop("The second argument must be a data frame")
  if (length(dim(ratetable))!=3) stop("The life table must have 3 dimensions in the (age, sex, birth year) order")
  if (dim(ratetable)[2]!=2) stop("The life table must have 3 dimensions in the (age, sex, birth year) order")
  if(!(dist %in% c("exponential","weibull","genweibull")))  stop("Argument 'dist' must be 'exponential', 'weibull', or 'genweibull' ")
  
  ####### data management
  
  time <- data[,as.character(formula[[2]][2])]
  event <- data[,as.character(formula[[2]][3])]
  cova <- as.matrix(data[,gsub("\\+", "", attr(terms(formula), "term.labels"))])
  age <- data[,age]
  sex <- data[,sex]
  year <- data[,year]
  if(!is.numeric(age)) stop("Argument 'age' must be a numeric vector")
  if(!is.character(sex)) stop("Argument 'sex' must be a character string")
  if(min(names(table(sex)) %in% c("female","male"))==0) stop("Argument 'p.sex' must be 'male' or 'female'")
  for(i in colnames(cova)){ if(!is.numeric(cova[,i])) stop("All covariates must be numeric")  }
  
  d <- cbind(time, event, cova, age, 1*(sex=="male"), year)
  na <- !is.na(apply(d, MARGIN=1, FUN = "sum"))
  
  time <- time[na]
  event <- event[na]
  cova <- cova[na,]
  age <- age[na]
  sex <- sex[na]
  year <- year[na]
  
  ###### Compute the expected mortality of individuals
  
   hP <- expectedhaz(ratetable, age, sex, year, time) # compute instantaneous hazards
   HP <- expectedcumhaz(ratetable, age, sex, year, time) # compute cumulative hazards
  
   ###### log likelihood functions
   
   logll1 <- function(sigma, nu, theta, beta, time, event, cova, hP, HP){
   return(-1*sum(event*log(hP+exp(cova%*%beta)*(
     (1/theta)*(1+(time/sigma)^nu)^((1/theta)-1)*(nu/sigma)*(time/sigma)^(nu-1) ) )
     - HP + exp(cova%*%beta)*(1-(1+(time/sigma)^nu)^(1/theta)))) }
   
   logll0 <- function(sigma, nu, theta, time, event, hP, HP){
     return(-1*sum(event*log(hP+(
       (1/theta)*(1+(time/sigma)^nu)^((1/theta)-1)*(nu/sigma)*(time/sigma)^(nu-1) ) )
       - HP + (1-(1+(time/sigma)^nu)^(1/theta)))) }
   
   if (dist=="genweibull"){
     label <- c("log sigma", "log nu", "log theta", gsub("\\+", "", attr(terms(formula), "term.labels")))
     
     init1 <- c(rep(0,3),rep(0,dim(cova)[2]))
    
     loglik1 <- function(par, time, event, cova, hP, HP){
          sigma <- exp(par[1]); nu <- exp(par[2]);  theta <- exp(par[3]); beta <- par[4:(3+dim(cova)[2])] 
          return(logll1(sigma, nu, theta, beta, time, event, cova, hP, HP)) }
     
     init0 <- rep(0,3)
     
     loglik0 <- function(par, time, event, cova, hP, HP){
       sigma <- exp(par[1]); nu <- exp(par[2]);  theta <- exp(par[3]) 
       return(logll0(sigma, nu, theta, time, event, hP, HP)) }
    }
  
   
   if (dist=="weibull"){
     label <- c("log sigma", "log nu", gsub("\\+", "", attr(terms(formula), "term.labels")))
     
     init1 <- c(rep(0,2),rep(0,dim(cova)[2]))
     
     loglik1 <- function(par, time, event, cova, hP, HP){
       sigma <- exp(par[1]); nu <- exp(par[2]);  theta <- 1; beta <- par[3:(2+dim(cova)[2])] 
       return(logll1(sigma, nu, theta, beta, time, event, cova, hP, HP)) }
     
     init0 <- rep(0,2)
     
     loglik0 <- function(par, time, event, cova, hP, HP){
       sigma <- exp(par[1]); nu <- exp(par[2]);  theta <- 1
       return(logll0(sigma, nu, theta, time, event, hP, HP)) }
   }
   
   if (dist=="exponential"){
     label <- c("log sigma", gsub("\\+", "", attr(terms(formula), "term.labels")))
     
     init1 <- c(rep(0,1),rep(0,dim(cova)[2]))
     
     loglik1 <- function(par, time, event, cova, hP, HP){
       sigma <- exp(par[1]); nu <- 1;  theta <- 1; beta <- par[2:(1+dim(cova)[2])] 
       return(logll1(sigma, nu, theta, beta, time, event, cova, hP, HP)) }
     
     init0 <- rep(0,1)
     
     loglik0 <- function(par, time, event, cova, hP, HP){
       sigma <- exp(par[1]); nu <- 1;  theta <- 1
       return(logll0(sigma, nu, theta, time, event, hP, HP)) }
   }
   
   
   logllmax0 <- optim(par = init0, fn = loglik0, time = time, event = event,
                      hP = hP, HP = HP )
   
   indic <- 0
   while(indic <= 5){
     ll_val <- logllmax0$value
     logllmax0 <- optim(par = logllmax0$par, fn = loglik0, time = time, event = event,
                        hP = hP, HP = HP )
     delta <- ll_val - logllmax0$value
     if(delta ==0) {indic = indic + 1}
   }
   
   logllmax0 <- optim(par = logllmax0$par, fn = loglik0, time = time, event = event,
                      cova = cova, hP = hP, HP = HP, hessian = TRUE)
   
   logllmax1 <- optim(par = init1, fn = loglik1, time = time, event = event,
                     cova = cova, hP = hP, HP = HP )
   
  indic <- 0
  while(indic <= 5){
    ll_val <- logllmax1$value
    logllmax1 <- optim(par = logllmax1$par, fn = loglik1, time = time, event = event,
                cova = cova, hP = hP, HP = HP )
    delta <- ll_val - logllmax1$value
    if(delta ==0) {indic = indic + 1}
  }
  
  logllmax1 <- optim(par = logllmax1$par, fn = loglik1, time = time, event = event,
                    cova = cova, hP = hP, HP = HP, hessian = TRUE)
  
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



