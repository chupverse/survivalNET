\name{splinecubeP}
\alias{splinecubeP}

\title{
The derivative of a Restricted Cubic Spline}
\description{
This helper function used in survivalFLEXNET allows to compute the values of the derivative 
of a restricted cubic spline.
}
\usage{
splinecubeP(time, gamma, m, mpos = NULL)
}

\arguments{
\item{time}{A vector of time where the spline is computed.}
\item{gamma}{A vector of the spline coefficients. The length of the coefficients vector
 should be equal to the total number of knots (m+2)}
\item{m}{A numeric value of the number of internal knots wanted. The first and last knots are 
automatically included, so there will be m+2 knots. The total number of knots should equal
the length of the coefficient vector.}
\item{mpos}{A vector of the different knots positions. The positions are the percentiles and must be
between ] 0:1[. The first and last knots positions will be automatically chosen to be at the 
0 and 1 percentile. Default is \code{NULL}, and if so, the knots positions will be chosen evenly through percentiles of log(times).}
}

\details{
The Restricted Cubic Spline is defined as :
\deqn{ s(x,\gamma) = \gamma_0 + \gamma_1*x + \sum_{j=2}^{m+1} \gamma_j*\nu_j}
with \deqn{\nu_j = (x-k_j)_+^3+\phi_j(x-k_{min})_+^3+(1-\phi_j)(x-k_{max})_+^3}
and \deqn{\phi_j = \frac{k_{max}-k_j}{k_{max}-k_{min}}}
where k are the knots positions and x = log(times).
\\
The derivative is defined as :
\deqn{ s'(x,\gamma) = \gamma_1 + \sum_{j=2}^{m+1} \gamma_j*\nu'_j}
with \deqn{\nu'_j = 3(x-k_j)_+^2+3*\phi_j(x-k_{min})_+^2+3*(1-\phi_j)(x-k_{max})_+^2}
and \deqn{\phi_j = \frac{k_{max}-k_j}{k_{max}-k_{min}}}

}

\value{
\item{spln}{The derivated spline values at the input times.}
\item{knots}{The position of the knots (including the first and last knot). Note that those are log(times), as the log-transformation is needed to compute the spline.}
\item{phi}{The \eqn{\phi_j} coefficients calculated with the knots positions.}
\item{nu}{The values of \eqn{\nu'_j} for each log(times).}
}

\examples{

times = c(1,10,20,40)
gammas = c(0.4,1.3,0.9,1.4)
m = 2

splP = splinecubeP(times, gammas, m, mpos = NULL)

splP$spln

}

\keyword{Spline}
