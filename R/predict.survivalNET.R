
predict.survivalNET <- function(object, type="relative", newdata=NULL, newtimes=NULL,
                                ratetable = NULL, ...){
  
if(is.null(newtimes))  { newtimes <- 0:max(object$y[,1]) }
  
if(!is.null(newdata))
  {
    if(!is.data.frame(newdata)) stop("Argument 'newdata' must be a data frame")
    indic <- gsub("\\+", "", attr(terms(object$formula), "term.labels") ) %in% names(newdata) 
    if( sum(!indic) > 0 ) stop("Missing predictor in the data frame")
    covariates <- as.matrix(newdata[, gsub("\\+", "", attr(terms(object$formula), "term.labels"))])
    
    if(type=="overall"){
      if(is.null(ratetable))stop("For overall survival, the 'ratetable' argument is necessary.")
      age = object$call$age
      year = object$call$year
      sex = object$call$sex
      indic <- c(age,year,sex, gsub("\\+", "", attr(terms(object$formula), 
                    "term.labels") ) ) %in% names(newdata)
      if( sum(!indic) > 0 ) stop("Missing predictor in the data frame.
                                 For overall suvival, newdata also needs
                                 'age', 'sex' and 'year' for the ratetable")
      covariates <- as.matrix(newdata[, gsub("\\+", "", attr(terms(object$formula), "term.labels"))])
      
    }
     }
  
if(is.null(newdata))  { 
  if(type=="overall"){
    
  if(is.null(ratetable))stop("For overall survival, the 'ratetable' argument is necessary.")
  age = object$asy$age
  year = object$asy$year
  sex = object$asy$sex
  }
  covariates <- object$x  }
  
  
if(object$dist=="genweibull")  {
    sigma <- exp(object$coefficients[1])
    nu <- exp(object$coefficients[2])
    theta <- exp(object$coefficients[3])
    beta <- object$coefficients[4:(3+dim(object$x)[2])]
    }
  
if(object$dist=="weibull")  {
    sigma <- exp(object$coefficients[1])
    nu <- exp(object$coefficients[2])
    theta <- 1
    beta <- object$coefficients[3:(2+dim(object$x)[2])]
    }
  
if(object$dist=="exponential")  {
    sigma <- exp(object$coefficients[1])
    nu <- 1
    theta <- 1
    beta <- object$coefficients[2:(1+dim(object$x)[2])]
    }
  
if(type=="relative") {
  fun <- function(x) { exp( exp(covariates%*%beta)*(1-(1+(x/sigma)^nu)^(1/theta)) ) }
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
predictions <- sapply(newtimes, FUN = "fun")

return(list(times=newtimes, predictions=predictions))
}

# predi <- predict(mod1, newdata = data.frame(age=c(50,60), sex01=c(0,1), stade=c(1,1)))
# plot(predi$times, predi$predictions[1,], col=2, xlab = "Time in years", ylab = "Net survival")
# points(predi$times, predi$predictions[2,], col=1)




