
predict.survivalFLEXNET <- function(object, type="survival", newdata=NULL, newtimes=NULL, ...){
  
if(is.null(newtimes))  { newtimes <- 0:max(object$y[,1]) }
  
if(!is.null(newdata))
  {
    if(!is.data.frame(newdata)) stop("Argument 'newdata' must be a data frame")
    indic <- gsub("\\+", "", attr(terms(object$formula), "term.labels") ) %in% names(newdata) 
    if( sum(!indic) > 0 ) stop("Missing predictor in the data frame")
    covariates <- as.matrix(newdata[, gsub("\\+", "", attr(terms(object$formula), "term.labels"))])
  }
  
if(is.null(newdata))  {  covariates <- object$x  }
  
      beta <- object$coefficients[1:dim(object$x)[2]]
      gamma <- object$coefficients[(1+dim(object$x)[2]):length(object$coefficient)]
        
if(type=="survival") {
  fun <- function(x) { exp(-exp(covariates%*%beta)+splinecube(time = x, 
                       gamma = gamma, m = object$m, mpos = object$mpos))  }
  }
  
if(type=="lp") {
  fun <- function(x) { covariates %*% beta }
  }
  
predictions <- sapply(newtimes, FUN = "fun")

return(list(times=newtimes, predictions=predictions))
}

# predi <- predict(mod1, newdata = data.frame(age=c(50,60), sex01=c(0,1), stade=c(1,1)))
# plot(predi$times, predi$predictions[1,], col=2, xlab = "Time in years", ylab = "Net survival")
# points(predi$times, predi$predictions[2,], col=1)




