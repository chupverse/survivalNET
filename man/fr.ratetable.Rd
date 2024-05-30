\name{fr.ratetable}
\alias{fr.ratetable}
\docType{data}
\title{
Expected Mortality Rates of the General French Population
}
\description{
An object of class ratetable for the expected mortality of the French population. It is an array with three dimensions: age, sex and year.
}
\usage{data(fr.ratetable)}
\format{
The format is "ratetable". The attributes are:
\tabular{rl}{
\code{dim} \tab A numeric vector with the length of each dimension. \cr											
\code{dimnames} \tab A character vector with the names of each variable of the three dimensions. \cr
\code{dimid} \tab A character vector with the identification of the dimensions: \code{age}, \code{year} and \code{sex}. \cr
\code{factor} \tab A numeric vector of indicators equals to 1 if the corresponding dimension does not vary \cr
\tab according to the time. Only the dimension related to \code{sex} equals 1. \cr
\code{cutpoints} \tab A list of the thresholds to identify the changes in the mortality rates according to \cr
\tab the time-dependent dimensions (NULL for \code{sex}). \cr
\code{class} \tab The class of the object: \code{ratetable}. \cr
}
}

\details{
The organization of a ratetable object is described in details by Therneau (1999) and Pohar (2006). The original data and updates can be downloaded from the Human Life-Table Database (HMD, The Human Mortality Database).
}

\source{
URL: www.mortality.org
}

\references{
Therneau and Offord. Expected Survival Based on Hazard Rates (Update), Technical Report, Section of Biostatistics, Mayo Clinic 63, 1999. 

Pohar and Stare. Relative survival analysis in R. Computer methods and programs in biomedicine, 81: 272-278, 2006.
}

\examples{
data(fr.ratetable)

is.ratetable(fr.ratetable)
}

\keyword{Ratetable}
