
predict.survivalNET <- function(object, type="relative", newdata=NULL, newtimes=NULL,
                                ratetable = NULL, ...){
  
if(is.null(newtimes))  { newtimes <- 0:max(object$y[,1]) }
  
if(!is.null(newdata))
  {
    if(!is.data.frame(newdata)) stop("Argument 'newdata' must be a data frame")
  
    covnames <- names(as.data.frame(sNET$x))
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
      ratetable_vars <- assign_ratetable_vars(unlist(lapply(ratetable_terms, extract_vars_within_function)))
      age = ratetable_vars$age
      year = ratetable_vars$year
      sex = ratetable_vars$sex
      indic <- c(age,year,sex, covnames)  %in% names(newdata) 
      if( sum(!indic) > 0 ) stop("Missing predictor in the data frame.
                                 For overall suvival, newdata also needs
                                 'age', 'sex' and 'year' for the ratetable")
      # covariates <- as.matrix(newdata[, covnames]) déja présent à la ligne 18 j'ai l'impression
      
    }
     }
  
if(is.null(newdata))  { 
  if(type=="overall"){
    
  if(is.null(ratetable))stop("For overall survival, the 'ratetable' argument is necessary.")
  age = object$ays$age
  year = object$ays$year
  sex = object$ays$sex
  }
  covariates <- object$x  ##rajouter strata
  }
  
  
  ### dist NET
  
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
  }
  
  ## m FLEXNET
  
  if("m" %in% names(object)){
    beta <- object$coefficients[1:(dim(object$x)[2])]
    gamma <- object$coefficients[(dim(object$x)[2]+1):(dim(object$t.table)[1])]
    m = object$m
    mpos = object$mpos
  }
  
  ### type 
  
if(type=="relative") {
  if ("dist" %in% names(object)) {
  fun <- function(x) { exp( exp(covariates%*%beta)*(1-(1+(x/sigma)^nu)^(1/theta)) ) }
  }
  if("m" %in% names(object)){
  fun <- function(x) { exp( -1*exp(splinecube(x, gamma, m, mpos)$spln)*exp(covariates%*%beta))}
  }
  
}
  
if(type=="lp") {
  fun <- function(x) { covariates %*% beta }
  }
  
if(type=="overall"){
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
  
  ###predictions 
  if ("dist" %in% names(object)) {
 predictions <- sapply(newtimes, FUN = "fun")
  }
  if ("m" %in% names(object)) {
 predictions <- fun(newtimes)
    }

return(list(times=newtimes, predictions=predictions))
}

# predi <- predict(mod1, newdata = data.frame(age=c(50,60), sex01=c(0,1), stade=c(1,1)))
# plot(predi$times, predi$predictions[1,], col=2, xlab = "Time in years", ylab = "Net survival")
# points(predi$times, predi$predictions[2,], col=1)




