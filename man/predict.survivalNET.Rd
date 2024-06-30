\name{predict.survivalNET}
\alias{predict.survivalNET}

\title{
Predicting Net Survival
}
\description{
Predict the net survival based on a model or algorithm from an object of the class \code{survivalNET}.
}

\usage{
\method{predict}{survivalNET}(object, type="survival", newdata, newtimes, ...)
}

\arguments{
  \item{object}{An object returned by the function \code{urvivalNET}.}
  \item{type}{An character sting indicating the type of predicted value: "lp" returns the linear predictor of the model related to the excess hazard, while "survival" (the default value) returns the net survival rates.}
  \item{newdata}{An optional data frame containing covariate values at which to produce predicted values. There must be a column for every covariate included in the training sample. The default value is \code{NULL}, the predicted values are computed for the subjects of the training sample.}
  \item{newtimes}{The times at which to produce predicted values. The default value is \code{NULL}, the predicted values are computed for the observed times in the training data frame.}
  \item{...}{For future methods.}
}

\value{
\item{times}{A vector of numeric values with the times of the \code{predictions}.}
\item{predictions}{A matrix with the predictions of survivals of each subject (lines) for each observed time (columns).}
}

\examples{
data(dataK) # the database with the observed sample
data(fr.ratetable) # the table with the expected mortality rates

# Estimate the multivariate net survival model 

sNET <- survivalNET(Surv(time, event) ~ biomarker + sex + stade + delay, data = dataK,
         ratetable=fr.ratetable, age="age", sex="sexchara", year="year", dist="weibull",
                    strata=NULL, weights=NULL)

# Consider women with a biomarker equals 0.5, a delayed diagnosis of an advanced cancer

sNET.pred <- predict(sNET, newdata = 
                       data.frame( biomarker=0.5, sex=2, stade=0, delay=1 ) )


# Plot the related net survival

plot(sNET.pred$times/365.24, sNET.pred$predictions, type = "l",
     xlab = "Time since the diagnosis in years", ylab="Predicted net survival")
}

\keyword{Predicting}
