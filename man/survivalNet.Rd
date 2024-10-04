\name{survivalNET}
\alias{survivalNET}

\title{
Parametric Net survival
}
\description{
This function allows to compute a parametric regression for modelling net survival.
}
\usage{
survivalNET(formula, data, ratetable, dist="weibull", weights=NULL)
}

\arguments{
\item{formula}{a formula object, with the response on the left of a ~ operator, and the terms on the right. The response must be a survival object as returned by the \code{Surv} function. On the right side, it is necessary to specify a '+\code{ratetable}(age = age, year = year, sex = sex)' term.}
\item{data}{A data frame in which to look for the variables included in the \code{formula}.}
\item{ratetable}{A table of expected rates of the studied event, such as \code{fr.ratetable} for the mortality of the French population, i.e. an array with three dimensions: age, sex and year.}
\item{dist}{A character string listed in the following table to define the baseline hazard function. By default, the Weibull distribution is used.}
\item{weights}{A numeric vector for a weighted likelihood, i.e. different contributions of individuals. By default, the weights equal 1.}
}

\details{
NOTE : For the ratetable, \code{age} at baseline must be in days, \code{year}, the baseline date in the study, must be of class 'date' under the date format: number of days since 01.01.1960 and \code{sex} must be of class 'character' with genders : 'male' or 'female'. 
\cr
\cr
It is possible to stratify on a categorical covariate by adding a '+\code{strata}(timevar)' term on the right side of the formula.
\cr
\cr
The following hazard functions available are:   
\tabular{rrl}{
exponential \tab Exponential distribution \tab \eqn{\lambda(t)=1/\sigma} \cr											
weibull \tab Weibull distribution \tab \eqn{\lambda(t)=\nu(\frac{1}{\sigma})^{\nu}t^{\nu-1}} \cr
gebweibull \tab Generalized Weibull distribution \tab \eqn{\lambda(t)=\frac{1}{\theta}\left(1+\left(\frac{t}{\sigma}\right)^{\nu}\right)^{\frac{1}{\theta}-1} \nu\left(\frac{1}{\sigma}\right)^{\nu} t^{\nu-1}}
}
with \eqn{\sigma}, \eqn{\nu},and \eqn{\theta>0}. 
}

\value{
\item{formula}{The \code{formula} object.}
\item{dist}{The distribution of the baseline hazard function of the excess hazard.}
\item{coefficients}{A vector with the regression coefficients of the excess hazard.}
\item{var}{The variance-covariance matrix of the regression \code{coefficients}.}
\item{t.table}{A data frame with the estimates, the standard deviations, Wald statistics, and p-values.}
\item{loglik}{A vector containing the log-likelihood of the model without and with covariates.}
\item{linear.predictors}{A vector with the linear predictor of the included individuals.}
\item{missing}{A vector of logical values indication whether individuals are included or excluded because of missing data.}
\item{n}{The number of individuals used in the fit.}
\item{nevent}{The number of events used in the fit.}
\item{y}{A matrix with two columns with the observed outcomes: \code{time} and \code{status}.}
\item{x}{A matrix with the observed covariates.}
\item{ays}{A data frame with age, birthyear of individuals and gender.}
\item{call}{Return the call with the specified arguments of the functon.}
}

\references{
Dickman PW, Sloggett A, Hills M, and Hakulinen T (2004). Regression models for relative survival. Stat Med. 23 (1):51-64. doi:10.1002/sim.1597.
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

sNET <- survivalNET(Surv(time, event) ~ stade + delay + ratetable(age = age, year = year,                            sex = sexchara), data = dataK, ratetable=fr.ratetable,
                         dist="weibull",weights=NULL)

sNET.pred <- predict(sNET, newdata=data.frame( stade=c(0,0,1,1), delay=c(0,1,0,1) ) )

for (i in 1:4) {lines(sNET.pred$times, sNET.pred$predictions[i,], col=2)}
}

\keyword{Regressing}
