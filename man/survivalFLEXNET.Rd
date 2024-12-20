\name{survivalFLEXNET}
\alias{survivalFLEXNET}

\title{
Flexible Parametric Net survival
}
\description{
This function allows to compute a flexible regression with splines for modelling net survival.
}
\usage{
survivalFLEXNET(formula, data, ratetable, m=3, mpos = NULL, init=NULL, delta_th=0,
weights=NULL)
}

\arguments{
\item{formula}{a formula object, with the response on the left of a ~ operator, and the terms on the
right. The response must be a survival object as returned by the \code{Surv} function. On the right
side, it is necessary to specify a '+\code{ratetable}(age = age, year = year, sex = sex)' term. The 
time must be in days.}
\item{data}{A data frame in which to look for the variables included in the \code{formula}.}
\item{ratetable}{A table of expected rates of the studied event, such as \code{fr.ratetable} for the
mortality of the French population, i.e. an array with three dimensions: age, sex and year.}
\item{m}{A numeric value of the number of internal knots wanted. The first and last knots are 
automatically included, so there will be m+2 knots.}
\item{mpos}{A vector of the different knots positions. The positions are the percentiles and must be
between \] 0:1\[. The first and last knots positions will be automatically chosen to be at the 
0 and 1 percentile. Default is NULL, and if so, the knots positions will be chosen evenly through percentiles of log(times).}
\item{init}{A numeric vector of the initial parameters values for the optimization process. Default
is set to NULL, which means that all parameters values will be set to 0.}
\item{delta_th}{A threshold value from which the optimization process stops. Default is set to 0. 
If 0, the difference between the log-likelihood values in the optimization process must be 0. 
If another value is chosen, the difference must be less than equal than this value. See 
\code{details} for the optimization process.}
\item{weights}{A numeric vector for a weighted likelihood, i.e. different contributions of 
individuals. By default, the weights equal 1.}
}

\details{
NOTE : For the ratetable, \code{age} at baseline must be in days, \code{year}, the baseline date in 
the study, must be of class 'date' under the date format: number of days since 01.01.1960 and 
\code{sex} must be of class 'character' with genders : 'male' or 'female'. 
\cr
\cr
It is possible to stratify on a categorical covariate by adding a '+\code{strata}(timevar)' term on
the right side of the formula.
\cr
\cr
The optimization process is a made so that, from initial values, the parameters are estimated with 
the \code{optim} function. The estimated parameters are then re-used as the initial parameters of
an optimization process, and while the difference between the old log-likelihood value and the new
log-likelihood value is different than 0 (or the value chosen with delta_th), the process continue.
When the difference value is below the threshold or equal to 0 five times, the process stops and the
optmization process is done one last time with those initial parameters to get the final estimated 
parameters.

}

\value{
\item{formula}{The \code{formula} object.}
\item{dist}{The distribution of the baseline hazard function of the excess hazard.}
\item{coefficients}{A vector with the regression coefficients of the excess hazard.}
\item{var}{The variance-covariance matrix of the regression \code{coefficients}.}
\item{t.table}{A data frame with the estimates, the standard deviations, Wald statistics, and 
p-values.}
\item{loglik}{A vector containing the log-likelihood of the model without and with covariates.}
\item{linear.predictors}{A vector with the linear predictor of the included individuals.}
\item{missing}{A vector of logical values indication whether individuals are included or excluded 
because of missing data.}
\item{n}{The number of individuals used in the fit.}
\item{nevent}{The number of events used in the fit.}
\item{y}{A matrix with two columns with the observed outcomes: \code{time} and \code{status}.}
\item{x}{A matrix with the observed covariates.}
\item{ays}{A data frame with age, birthyear of individuals and gender.}
\item{call}{Return the call with the specified arguments of the functon.}
}

\references{
Dickman PW, Sloggett A, Hills M, and Hakulinen T (2004). Regression models for relative survival. 
Stat Med. 23 (1):51-64. doi:10.1002/sim.1597.
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

sNET <- survivalNET(Surv(time, event) ~ stade + delay + ratetable(age = age, year = year, 
                    sex = sexchara), data = dataK, ratetable=fr.ratetable,
                    dist="weibull", init = NULL, delta_th = 0, weights = 0)

sNET.pred <- predict(sNET, newdata=data.frame( stade=c(0,0,1,1), delay=c(0,1,0,1) ) )

for (i in 1:4) {lines(sNET.pred$times, sNET.pred$predictions[i,], col=2)}
}

\keyword{Regressing}
