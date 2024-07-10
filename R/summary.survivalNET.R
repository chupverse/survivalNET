
summary.survivalNET <- function (object, digits = 4, conf.int = 0.95, scale = 1
                                 , ...)
{
  
  list <- c(call = object$formula, 
            n = object$n,
            loglik = c(object$loglik)
              )
  return(list)
  } 


