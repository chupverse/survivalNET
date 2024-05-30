\name{expectedhaz}
\alias{expectedhaz}

\title{
Expected Instantaneous Hazard
}
\description{
This function returns the expected instantaneous hazard from a ratetable according to age, sex and year at the baseline.
}
\usage{
expectedhaz(ratetable, age, sex, year, time)
}

\arguments{
\item{ratetable}{A table of expected rates of the studied event, such as \code{fr.ratetable} for the motality of the French population, i.e. an array with three dimensions: age, sex and year.}
\item{age}{The patient age in days at baseline.}
\item{sex}{A character vector with the gender: \code{male} or \code{female}.}
\item{year}{A numeric vector with the calendar year at baseline.}
\item{time}{A scalar with the time post-baseline in days.}
}

\examples{
data(fr.ratetable)

# the expected instantaneous hazard of death at 1-year post-inclusion
# for a French 50 years-old man included in a study in 2010
expectedhaz(fr.ratetable, age=50, sex="male", year=2010, time=365)
}

\keyword{Ratetable}
