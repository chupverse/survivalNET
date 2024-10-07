\name{plot.survivalNET}
\alias{plot.survivalNET}

\title{
Calibration Plot for a Relative Survival Model with Additive Hazards
}
\description{
A calibration plot of an object of the class \code{survivalNET}.
}

\usage{
\method{plot}{survivalNET}(x, n.groups=5, pro.time=NULL,
newdata=NULL, ratetable, age, year, sex, ...)
}

\arguments{
  \item{x}{An object returned by the function \code{survivalNET}.}
   \item{n.groups}{A numeric value with the number of groups by their class probabilities. The default is 5.}
  \item{pro.time}{The prognostic time at which the calibration plot of the survival probabilities.}
    \item{newdata}{An optional data frame containing the new sample for validation with covariate values, follow-up times, and event status. The default value is \code{NULL}: the calibration plot is performed from the same subjects of the training sample.}
\item{ratetable}{A table of expected rates of the studied event, such as \code{fr.ratetable} for the mortality of the French population, i.e. an array with three dimensions: age, sex and year.}
\item{age}{A numeric vector with the patient age at baseline in days.}
\item{year}{A numeric vector with the baseline date in the study respecting the date format: number of days since 01.01.1960.}
\item{sex}{A character vector with the gender: "male" or "female".}
  \item{...}{Additional arguments affecting the plot.}
}

\details{
The plot represents the non-parametric net survival obtained by the Pohar-Perme estimator and the related 95\% confidence intervals, against the mean of the predictive values for individuals stratified into groups of the same size according to the percentiles. The identity line is usually included for reference.
}

\value{No return value for this S3 method.}

\examples{
data(dataK) # the database with the observed sample
data(fr.ratetable) # the table with the expected mortality rates

# The multivariate PH regression with a Weibull distribution

sNET <- survivalNET(Surv(time, event) ~ biomarker + sex + stade + delay +
                    ratetable(age = age, year = year, sex = sexchara),
                    data = dataK, ratetable=fr.ratetable,
                    dist="weibull")
                    
# the calibration plots for predictions at 6

plot(sNET, n.groups=3, pro.time=6*365.24, 
  ratetable=fr.ratetable, age="age", sex="sexchara", year="diagnum",
  main="Predction at 6 years")
}

\seealso{
\code{\link{plot.default}}, \code{\link{survivalNET}}.
}

\keyword{Plotting}
