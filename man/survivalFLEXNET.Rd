\name{survivalFLEXNET}
\alias{survivalFLEXNET}

\title{
Flexible Parametric Net survival
}
\description{
This function allows to compute a flexible regression with splines for modelling net survival.
}
\usage{
survivalFLEXNET(formula, data, ratetable, m=3, mpos = NULL, mquant = NULL, init=NULL, delta_th=0,
weights=NULL, m_s = NULL, mpos_s = NULL, mquant_s = NULL, Kref = NULL)
}

\arguments{
\item{formula}{a formula object, with the response on the left of a ~ operator, and the terms on the
right. The response must be a survival object as returned by the \code{Surv} function. On the right
side, it is necessary to specify a '+\code{ratetable}(age = age, year = year, sex = sex)' term. The 
time must be in days.}
\item{data}{A data frame in which to look for the variables included in the \code{formula}.}
\item{ratetable}{A table of expected rates of the studied event, such as \code{fr.ratetable} for the
mortality of the French population, i.e. an array with three dimensions: age, sex and year.}
\item{m}{A numeric value of the number of internal knots wanted. With the boundary knots, there will
be m+2 knots.}
\item{mpos}{A vector of the different m+2 knots values positions. The values must be specified on 
the log scale, as no transformation is done on those values compared to the time vector provided. 
Default is \code{NULL}, and if so, the knots positions will be chosen evenly
through percentiles of log(times).}
\item{mquant}{A vector of the different m+2 knots quantile positions. The positions are the 
percentiles and must be between [0:1]. Default is \code{NULL}, and if so, the knots positions will 
be chosen evenly through percentiles of log(times). Will be overwritten by \code{mpos} if both are
specified.}
\item{init}{A numeric vector of the initial parameters values for the optimization process. Default
is set to NULL, which means that all parameters values will be set to 0.}
\item{delta_th}{A threshold value from which the optimization process stops. Default is set to 0. 
If 0, the difference between the log-likelihood values in the optimization process must be 0. 
If another value is chosen, the difference must be less than equal than this value. See 
\code{details} for the optimization process.}
\item{weights}{A numeric vector for a weighted likelihood, i.e. different contributions of 
individuals. By default, the weights equal 1.}
\item{m_s}{A numeric value of the number of internal knots wanted for the splines for each level of stratification. With the boundary knots, there will
be m_s+2 knots.}
\item{mpos_s}{A vector of the different m_s+2 knots values positions for the splines for each level of stratification. The values must be specified on 
the log scale, as no transformation is done on those values compared to the time vector provided. 
Default is \code{NULL}, and if so, the knots positions will be chosen evenly
through percentiles of log(times).}
\item{mquant_s}{A vector of the different m_s+2 knots quantile positions for the splines for each level of stratification. The positions are the 
percentiles and must be between [0:1]. Default is \code{NULL}, and if so, the knots positions will 
be chosen evenly through percentiles of log(times). Will be overwritten by \code{mpos_s} if both are
specified.}
\item{Kref}{The level of reference when stratifying on a covariate. Default value is NULL. If NULL the value taken will be the maximum one (if numeric) or the last alphabetically (if character string).}
}

\details{
NOTE : For the ratetable, \code{age} at baseline must be in days, \code{year}, the baseline date in 
the study, can be of class 'numeric', 'Date', 'POSIXt' or 'date' under the date format: number of days since 01.01.1960 and 
\code{sex} must be of class 'character' with genders : 'male' or 'female'. 
\cr
\cr
It is possible to stratify on a categorical covariate by adding a '+\code{strata}(timevar)' term on
the right side of the formula.
\cr
\cr
The optimization process is a made so that, from initial values, the parameters are estimated with 
the \code{optim} function. The estimated parameters are then re-used as the initial parameters of
an optimization process, and while the difference between the old log-likelihood value and the new
log-likelihood value is different than 0 (or the value chosen with delta_th), the process continue.
When the difference value is below the threshold or equal to 0 five times, the process stops and the
optmization process is done one last time with those initial parameters to get the final estimated 
parameters.

}

\value{
an object of class "\code{survivalNET}".
\item{formula}{The \code{formula} object.}
\item{coefficients}{A vector with the regression coefficients of the excess hazard.}
\item{var}{The variance-covariance matrix of the regression \code{coefficients}.}
\item{t.table}{A data frame with the estimates, the standard deviations, Wald statistics, and 
p-values.}
\item{loglik}{A vector containing the log-likelihood of the model without and with covariates.}
\item{linear.predictors}{A vector with the linear predictor of the included individuals.}
\item{missing}{A vector of logical values indication whether individuals are included or excluded 
because of missing data.}
\item{n}{The number of individuals used in the fit.}
\item{nevent}{The number of events used in the fit.}
\item{y}{A matrix with two columns with the observed outcomes: \code{time} and \code{status}.}
\item{x}{A matrix with the observed covariates.}
\item{ays}{A data frame with age, birthyear of individuals and gender.}
\item{m}{A numeric value of the number of internal knots specified by the user.}
\item{mpos}{The position of the knots, either provided by the user or computed automatically.}
\item{mquant}{The quantile values used for the knots selection, either provided by the user or
computed automatically.}
\item{correstab}{Correspondance table of the stratified term values.}
\item{xlevels}{The different levels of the stratified term.}
\item{levelsval}{The stratified term individual values.}

}

\references{
Royston P, Parmar MK. Flexible parametric proportional-hazards and proportional-odds models for censored survival data, with application to prognostic modelling and estimation of treatment effects. Stat Med. 2002 Aug 15;21(15):2175-97. doi: 10.1002/sim.1203. PMID: 12210632.
}

\examples{
data(dataK) # the database with the observed sample
data(fr.ratetable) # the table with the expected mortality rates

# The non-parametric estimation of the net survival curves 

PPest <- rs.surv(Surv(time, event) ~ stade + delay, data = dataK,
                   ratetable = fr.ratetable, method = "pohar-perme",
                   rmap = list(age = age, sex = sex, year = year))

plot(PPest) 

# The flexible parametric estimation of the net survival curves 

sFLEX <- survivalFLEXNET(Surv(time, event) ~ stade + delay + ratetable(age = age, year = year, 
                    sex = sexchara), data = dataK, ratetable=fr.ratetable,
                    m = 2, init = NULL, delta_th = 0, weights = NULL, m_s = NULL, mpos_s = NULL, mquant_s = NULL, Kref = NULL)

sFLEX.pred <- predict(sFLEX, newdata=data.frame( stade=c(0,0,1,1), delay=c(0,1,0,1) ) )

for (i in 1:4) {lines(sFLEX.pred$times, sFLEX.pred$predictions[i,], col=2)}
}

\keyword{Regressing}
