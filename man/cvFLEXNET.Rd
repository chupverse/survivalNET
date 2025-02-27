\name{cvFLEXNET}
\alias{cvFLEXNET}

\title{
Cross-validation method for a Flexible Relative Survival Model with Additive Hazards
}
\description{
Cross-validation method for the hyper-parameters of an object of class \code{survivalFLEXNET}.
}

\usage{
cvFLEXNET(formula, pro.time=NULL, data, ratetable, cv=10, 
                             m = 2, mpos = NULL, init = NULL,
                             delta_th = 0, weights = NULL)
}

\arguments{
  \item{formula}{The \code{formula} object.}
  \item{pro.time}{The prognostic time at which the metric is evaluated. Default value is \code{NULL} and if so, the median 
  of the times in chosen.}
  \item{data}{A data frame in which to look for the variables included in the \code{formula}.}
  \item{ratetable}{A table of expected rates of the studied event, such as \code{fr.ratetable} for
  the mortality of the French population, i.e. an array with three dimensions: age, sex and year.}
  \item{cv}{The number of splits for cross-validation. The default value is 10.}
  \item{m}{A numeric value of the number of internal knots wanted. The first and last knots are 
automatically included, so there will be m+2 knots.}
  \item{mpos}{A vector of the different knots positions. The positions are the percentiles and must be
between \] 0:1\[. The first and last knots positions will be automatically chosen to be at the 
0 and 1 percentile. Default is \code{NULL}, and if so, the knots positions will be chosen evenly through percentiles of log(times).}
\item{init}{A numeric vector of the initial parameters values for the optimization process. Default
is set to NULL, which means that all parameters values will be set to 0.}
\item{delta_th}{A threshold value from which the optimization process stops. Default is set to 0. 
If 0, the difference between the log-likelihood values in the optimization process must be 0. 
If another value is chosen, the difference must be less than equal than this value. See \code{details} of \code{survivalFLEXNET} for the optimization process.}
\item{weights}{A numeric vector for a weighted likelihood, i.e. different contributions of 
individuals. By default, the weights equal 1.}
}

\details{The hyper-parameters are "m" and "mpos". It is only possible to do cross-valiation on one at the time. The user must input either multiple 'm' with 'mpos = NULL' or only one 'm' and multiple differents knots positions."}

\value{
\item{optimal}{A list giving the optimal value for each parameter.}
\item{results}{A data frame listing the parameters combinaison and their metrics values.}
}

\examples{
data(dataK) # the database with the observed sample
data(fr.ratetable) # the table with the expected mortality rates

tune.FLEX = cvFLEXNET(Surv(time, event) ~ stade + delay + 
    ratetable(age=age, year=year, sex=sexchara), data = dataK[1:99,],
    ratetable=fr.ratetable, cv = 2, m = c(1,2))

tune.FLEX$optimal

tune.FLEX$results 
}

\seealso{
\code{\link{survivalFLEXNET}}.
}

\keyword{Cross-validation}
