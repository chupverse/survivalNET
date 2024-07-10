\name{predict.survivalNET}
\alias{predict.survivalNET}

\title{
Predicting Net Survival
}
\description{
Predict the net survival based on a model or an algorithm from an object of the class \code{survivalNET}.
}

\usage{
\method{predict}{survivalNET}(object, type="survival", newdata, newtimes, ...)
}

\arguments{
  \item{object}{An object returned by the function \code{urvivalNET}.}
  \item{type}{A character string indicating the type of predicted value : "lp" returns the linear predictor of the model related to the excess hazard, while "survival" (the default value) returns the net survival rates.}
  \item{newdata}{An optional data frame containing covariate values at which to produce predicted values. There must be a column for every covariate included in the training sample. If \code{NULL} (the default value), the predicted values are computed for the subjects of the training sample.}
  \item{newtimes}{The times at which to produce predicted values. If \code{NULL} (the default value), the predicted values are computed for the observed times in the training data frame.}
  \item{...}{For future methods.}
}

\value{
\item{times}{A vector of numeric values with the times of the \code{predictions}.}
\item{predictions}{A matrix with the predictions of survivals of each subject (lines) for each observed time (columns).}
}

\examples{
data(dataK) # the database with the observed sample
data(fr.ratetable) # the table with the expected mortality rates

# The non-parametric estimation of the net survival curves 

PPest <- rs.surv(Surv(time, event) ~ stade + delay, data = dataK,
                   ratetable = fr.ratetable, method = "pohar-perme",
                   rmap = list(age = age, sex = sex, year = year))

plot(PPest) 

# The parametric estimation of the net survival curves (Weibull distribution)

sNET <- survivalNET(Surv(time, event) ~ stade + delay, data = dataK, ratetable=fr.ratetable,
                     age="age", sex="sexchara", year="year", dist="weibull",
                     strata=NULL, weights=NULL)

sNET.pred <- predict(sNET, newdata=data.frame( stade=c(0,0,1,1), delay=c(0,1,0,1) ) )

for (i in 1:4) {lines(sNET.pred$times, sNET.pred$predictions[i,], col=2)}
}

\keyword{Predicting}
