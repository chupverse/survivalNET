\name{dataK}
\alias{dataK}
\docType{data}
\title{
Observed Mortality Data of Patients with Cancer
}
\description{
An data frame with a simulated sample of patients with a cancer and their follow-up.
}
\usage{data(dataK)}
\format{
The format is "data.frame". The names of the columns are:
\tabular{rl}{
"sex"       "age"       "year"      "stade"     "delay"     "biomarker" "sexchara"  "event"     "time"    


\code{sex} \tab A numeric vector equals 1 for male and 2 for male. \cr											
\code{age} \tab A  numeric vector with the age of the patient at the time of cancer diagnosis (baseline). \cr
\code{year} \tab A numeric vector with the date of diagnosis (in date format, i.e. the number of days since 1Jan60). \cr
\code{stade} \tab A numeric vector with the disease stage: 1 for advanced cancer and 0 otherwise. \cr
\tab according to the time. Only the dimension related to \code{sex} equals 1. \cr
\code{delay} \tab A numeric vector indicating if the cancer was early diagnosed: 1 for delayed diagnosis and 0 otherwise. \cr
\code{biomarker} \tab  A numeric vector with the level of a continuous biomarker with is associated with an increase cancer-specific mortality. \cr
\code{sexchara} \tab  A character vector with the patient gender: "male" or "female". \cr
\code{event} \tab  A numeric vector equals 1 for death and 0 for censoring. \cr
\code{time} \tab  A numeric vector with the follow-up time in days since the cancer diagnostic.
}
}

\details{
The data frame was obtained by simulations. The French mortality tables were used for the expected mortality and a proportional hazard model with an Exponential distribution for the baseline hazard.
}

\examples{
data(dataK)

# Plotting the observed survival by using the Kaplan-Meier estimator

plot(survfit(Surv(time/365.241, event) ~ 1, data = dataK),
     ylab="Patient survival", xlab="Post-diagnostic time in years")
}

\keyword{Data}