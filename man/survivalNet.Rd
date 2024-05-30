\name{survivalNET}
\alias{survivalNET}

\title{
Parametric Net survival
}
\description{
This function allows to compute a parametric regression for modelling net survival.
}
\usage{
survivalNET(formula, data, ratetable, age, sex, year, dist="weibull", strata=NULL, weights=NULL)
}

\arguments{
\item{formula}{a formula object, with the response on the left of a ~ operator, and the terms on the right. The response must be a survival object as returned by the \code{Surv} function.}
\item{data}{A data frame in which to look for the variables included in the \code{formula}.}
\item{ratetable}{A table of expected rates of the studied event, such as \code{fr.ratetable} for the motality of the French population, i.e. an array with three dimensions: age, sex and year.}
\item{age}{A numeric vector with the patient ages in days at baseline.}
\item{sex}{A character vector with the genders: \code{male} or \code{female}.}
\item{year}{A numeric vector with the entry dates in the study respecting the date format, i.e. in number of days since 01.01.1960.}
\item{dist}{One of the character strings listed in the following table to define the baseline hazard function. By default, the Weibull distribution is used.}
\item{strata}{An optional argument indicating the strata for which separate baseline hazard functions are fit.}
\item{weights}{A numeric vector with a weigthed likelihood, i.e. different contributions of individuals. By defaut, the weigths equal 1.}
}

\details{
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
\item{loglik}{A vector of length 2 containing the log-likelihood of the model without and with the ccovariates.}
\item{linear.predictors}{A vector with the linear predictor of the included individuals.}
\item{missing}{A vector of logical values indication whether individuals are included or excluded because of missing data.}
\item{n}{The number of individuals used in the fit.}
\item{nevent}{The number of events used in the fit.}
\item{y}{A matrix with two columns with the observed outcomes: \code{time} and \code{status}.}
\item{x}{A matrix with the observed covariates.}
\item{asy}{A data frame with age, gender, and birthyear of individuals}
}

\references{
Dickman PW, Sloggett A, Hills M, and Hakulinen T (2004). Regression models for relative survival. Stat Med. 23 (1):51-64. doi:10.1002/sim.1597.
}

\examples{
####### Observed survival

Sfun <- function(time, sigma, nu, theta, covariate, coef, ratetable, age, sex, year)
  {
  exp(-expectedcumhaz(ratetable, age, sex, year, time)) *
  exp((1-(1+(time/sigma)^nu)^(1/theta)) * exp(covariate \%*\% coef))
  }

######## function for times-to-event simulation

SimT <- function(n, logsigma, lognu, logtheta, covariate, coef, ratetable, age, sex, year)
  {
  survivaltimes = sapply(1:n,
                    function(i) {
                     u <- runif(1)
                     time <- uniroot(function(x) 1- Sfun(x, exp(logsigma), exp(lognu), exp(logtheta),
                              covariate[i,], coef, ratetable,
                              age[i], sex[i], year[i]) - u, interval = c(0, 1000))$root
                   return(time)
                                  }
                         )
  return(survivaltimes)
  }

######## function times-to-censoring simulation

SimC <- function(u, sigma) { (1/-sigma) * log(1-u) }

######## Simulations

n <- 2000 #sample size

# Patient characteristics at baseline
sex01 <- rbinom(n, 1, 0.5)
sex <- rep("male", n)
sex[sex01==1] <- "female"
age <- floor(runif(n, 0, 78)) 
birthyear <- floor(runif(n, 1940, 1981))
stade <- rbinom(n, 1, 0.3)
bruit <- rnorm(n, 0, 1)

data("fr.ratetable")

p <- c(1.25, 0.9, 0, 0.6, 0.4, 0.01)
timeT <- SimT(n, logsigma = p[1], lognu = p[2], logtheta = p[3], covariate = cbind(sex01, stade, age),
                 coef = c(p[4], p[5], p[6]), ratetable = fr.ratetable,
                 age = age, sex = sex, year = birthyear)

timeC <- SimC(runif(n), 0.2) #censoring time

simdata <- data.frame(
  time = pmin(timeT, timeC),
  status = 1*(timeC >= timeT),
  sex, sex01, stade, age, birthyear)

mod1 <- survivalNET(Surv(time, status) ~ sex01 + stade + age, data = simdata, ratetable=fr.ratetable,
              age="age", sex="sex", year="birthyear", dist="weibull")
}

\keyword{Net survival}
