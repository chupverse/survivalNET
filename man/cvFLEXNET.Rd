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
                             m = 2, mpos = NULL, mquant = NULL, init = NULL,
                             delta_th = 0, weights = NULL, m_s = NULL, Kref = NULL, metric = "ibs")
}

\arguments{
\item{formula}{The \code{formula} object.}
\item{pro.time}{The prognostic time at which the metric is evaluated. Default value is \code{NULL} and if so, the median 
of the times in chosen.}
\item{data}{A data frame in which to look for the variables included in the \code{formula}.}
\item{ratetable}{A table of expected rates of the studied event, such as \code{fr.ratetable} for
the mortality of the French population, i.e. an array with three dimensions: age, sex and year.}
\item{cv}{The number of splits for cross-validation. The default value is 10.}
\item{m}{A list of numeric values of the number of internal knots wanted. With the boundary knots, there will
be m+2 knots.}
\item{mpos}{A list of vectors of the different m+2 knots values positions. The values must be specified on 
the log scale, as no transformation is done on those values compared to the time vector provided. 
Default is \code{NULL}, and if so, the knots positions will be chosen evenly
through percentiles of log(times).}
\item{mquant}{A list of vectors of the different m+2 knots quantile positions. The positions are the 
percentiles and must be between [0:1]. Default is \code{NULL}, and if so, the knots positions will 
be chosen evenly through percentiles of log(times). Will be overwritten by \code{mpos} if both are
specified.}
\item{init}{A numeric vector of the initial parameters values for the optimization process. Default
is set to NULL, which means that all parameters values will be set to 0.}
\item{delta_th}{A threshold value from which the optimization process stops. Default is set to 0. 
If 0, the difference between the log-likelihood values in the optimization process must be 0. 
If another value is chosen, the difference must be less than equal than this value. See \code{details} of \code{survivalFLEXNET} for the optimization process.}
\item{weights}{A numeric vector for a weighted likelihood, i.e. different contributions of 
individuals. By default, the weights equal 1.}
\item{m_s}{A list of numeric values of the number of internal knots wanted for the splines for each level of stratification. With the boundary knots, there will
be m_s+2 knots.}
\item{Kref}{The level of reference when stratifying on a covariate. Default value is NULL. If NULL the value taken will be the maximum one (if numeric) or the last alphabetically (if character string).}
\item{metric}{The metric used to select the best set of parameters from the cross-validation. See the function 'metrics' for more details. Default is set to "ibs".}
}

\details{The hyper-parameters are 'm', 'mpos', 'mquant' and 'm_s'. For 'm', 'mpos' and 'mquant' it is only possible to do cross-valiation on one at the time. The user must input either multiple 'm' with 'mpos = NULL' and 'mquant = NULL' or only one 'm' and multiple differents knots 
positions or quantile values."}

\value{
\item{optimal}{A list giving the optimal value for each parameter.}
\item{results}{A data frame listing the parameters combinaison and their metrics values.}
}

\examples{
data(dataK) # the database with the observed sample
data(fr.ratetable) # the table with the expected mortality rates

tune.FLEX = cvFLEXNET(Surv(time, event) ~ stade + delay + 
    ratetable(age=age, year=year, sex=sexchara), data = dataK[1:99,],
    ratetable=fr.ratetable, cv = 2, m = c(1,2), mquant = NULL, init = NULL,
                             delta_th = 0, weights = NULL, m_s = NULL, Kref = NULL, metric = "ibs")

tune.FLEX$optimal

tune.FLEX$results 
}

\seealso{
\code{\link{survivalFLEXNET}}.
}

\keyword{Cross-validation}
