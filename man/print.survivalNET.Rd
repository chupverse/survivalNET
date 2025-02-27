\name{print.survivalNET}
\alias{print.survivalNET}

\title{
S3 Method for Printing a 'survivalNET' object
}
\description{
Print the main parameters of the estimated model.
}

\usage{
\method{print}{survivalNET}(x, digits=4, ...)
}

\arguments{
  \item{x}{An object returned by the function \code{survivalNET}.}
  \item{digits}{An optional integer for the number of digits to print when printing numeric values.}
  \item{...}{For future methods.}
}

\value{No return value for this S3 method.}

\seealso{
  \code{\link{survivalNET}}.
}

\examples{
data(dataK) # the database with the observed sample
data(fr.ratetable) # the table with the expected mortality rates

# The parametric estimation of the net survival curves (Weibull distribution)
sNET <- survivalNET(Surv(time, event) ~ stade + delay + ratetable(age=age, sex=sexchara, year=year),
data = dataK, ratetable=fr.ratetable, dist="weibull", weights=NULL)

print(sNET)
}

\keyword{Printing}
