\name{metrics}
\alias{metrics}

\title{
Metrics to Evaluate the Prognostic Capacities
}

\description{
Compute several metrics to evaluate the prognostic capacities with time-to-event data.
}

\usage{
metrics(formula, data, prediction.matrix, prediction.times, metric,
pro.time=NULL, ROC.precision=seq(.01, .99, by=.01))
}


\arguments{
  \item{formula}{A formula object, with a survival object \code{Surv(time,event)} on the left of a ~
  operator, and only "\code{~1}" on the right side.}
  \item{data}{A data frame for in which to look for the variables related to the status of the 
  follow-up time (\code{times}) and the event (\code{failures}).}
  \item{prediction.matrix}{A matrix or data.frame object of numeric values of dimension n*p, with n 
  being the number of suject and p being the number of prognostic times.}
  \item{prediction.times}{A vector of numeric values with the times of the \code{predictions} (same
  length than the number of columns of \code{prediction.matrix}).}
  \item{metric}{The metric to compute. See details.}
  \item{pro.time}{This optional value of prognostic time represents the maximum delay for which the 
  capacity of the variable is evaluated. The same unit than the one used in the argument times. Not
  used for the following metrics: "loglik", "ibs", "bll", and "ibll". Default value is the time at 
  which half of the subjects are still at risk.}
   \item{ROC.precision}{An optional argument with the percentiles (between 0 and 1) of the 
   prognostic variable used for computing each point of the time dependent ROC curve. Only used when
   \code{metric="auc"}. 0 (min) and 1 (max) are not allowed. By default, the precision is 
   \code{seq(.01,.99,.01)}.}
 }

\details{
The following metrics can be used: "bs" for the Brier score at the prognostic time \code{pro.time}, 
"ci" for the concordance index at the prognostic time \code{pro.time}, "loglik" for the 
log-likelihood, "ibs" for the integrated Brier score up to the last observed time of event, "ibll" 
for the integrated binomial log-likelihood up to the last observed time of event, "bll" for the 
binomial log-likelihood, "ribs" for the restricted integrated Brier score up to the prognostic time
\code{pro.time}, "ribll" for the restricted integrated binomial log-likelihood up to the last 
observed time of event, "bll" for the binomial log-likelihood, "auc" for the area under the 
time-dependent ROC curve up to the prognostic time \code{pro.time}.
}

\value{
A numeric value with the metric estimation.
}

\examples{

data(dataK)
data(fr.ratetable)
  
# The estimation of the model
model <- survivalFLEXNET(formula = Surv(time, event)~delay + stade + 
ratetable(age, year, sexchara), data = dataK, ratetable = fr.ratetable, m = 2)

#Predictions of the model

pred <- predict(model, newtimes = model$y[,1][model$y[,2] == 1]) 

# The apparent AUC at 10-year post-transplantation
metrics(formula = Surv(time, event)~1, data=dataK,
   prediction.matrix = pred$predictions, prediction.times=pred$times, metric
  ="auc", pro.time=10)

# The integrated Brier score up to 10 years post-transplanation
metrics(formula = Surv(time, event)~1, data=dataK,
  prediction.matrix = pred$predictions, prediction.times=pred$times, metric 
  ="ribs", pro.time=10)
}

\keyword{Metrics}
