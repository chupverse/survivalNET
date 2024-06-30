\name{expectedcumhaz}
\alias{expectedcumhaz}

\title{
Expected Cummulative Hazard
}
\description{
This function returns the expected cumulative hazard from a ratetable according to age, sex and year at the baseline.
}
\usage{
expectedcumhaz(ratetable, age, year, sex, time, method="exact", subdivisions = 100)
}

\arguments{
\item{ratetable}{A table of expected rates of the event, such as \code{fr.ratetable} for the motality of the French population, i.e. an array with three dimensions: age, sex and year. The names of the arguments \code{age} and \code{year} must be a sequence of years by steps equal to 1 years.}
\item{age}{A scalar with the patient age in days at baseline.}
\item{year}{A scalar with the baseline date  in days since 01.01.1960.}
\item{sex}{A character vector with the gender: "male" or "female".}
\item{time}{A scalar with the post-baseline time in days.}
\item{method}{The type of estimation: "exact" for the sum of the hazard rates per day during the periode from 0 to \code{time}, or "trapezoidal" for a numeric integratation using the composite trapezoidal rule.}
\item{subdivisions}{A scalar with the number of tapazes if \code{method="trapezoidal"}.}
}

\examples{
data(fr.ratetable)

# the expected cummulative hazard of death at 1-year post-inclusion
# for a French 50 years-old man included in a study in 2010

date.ddmmmyy(14000)
expectedcumhaz(fr.ratetable, age=50, sex="male", year=14000, time=365)
}

\keyword{Ratetable}
