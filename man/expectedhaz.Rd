\name{expectedhaz}
\alias{expectedhaz}

\title{
Expected Instantaneous Hazard
}
\description{
This function returns the expected instantaneous hazard from a rate table according to age, sex and year at baseline.
}
\usage{
expectedhaz(ratetable, age, year, sex, time)
}

\arguments{
\item{ratetable}{A table of expected rates of the event, such as \code{fr.ratetable} for the mortality of the French population, i.e. an array with three dimensions: age, sex and year. The names of the arguments \code{age} and \code{year} must be a sequence of years by steps equal to 1 years.}
\item{age}{A scalar with the patient age in days at baseline.}
\item{year}{A scalar with the baseline date  in days since 01.01.1960.}
\item{sex}{A character vector with the gender: "male" or "female".}
\item{time}{A scalar with the post-baseline time in days.}
}

\examples{
data(fr.ratetable)

# check that the steps equal to 1 years

attributes(fr.ratetable)$dimnames$age
attributes(fr.ratetable)$dimnames$year

# the expected instantaneous hazard of death at 1-year post-inclusion
# for a French 50 years-old man included in a study the 02.03.2005

date.ddmmmyy(14000)
expectedhaz(fr.ratetable, age=50*365, sex="male", year=2010, time=14000)
}
