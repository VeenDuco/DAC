\name{DAC.uniform}
\alias{DAC.uniform}
\title{DAC For Uniform Benchmark Prior And Expert Priors}
\usage{
DAC.uniform(from, to, by, data, priors, lb.bench, ub.bench, n.iter = 10000)
}
\description{
Calculates DAC values for multiple experts that provide prections in the form of a prior. The benchmark is set to be a uniform distribution with user specified boundaries.
}

\arguments{
\item{from}{Lower bound of the parameter space that is to be evaluated, as in the seq function of the base package.}
\item{to}{Upper bound of the parameter space that is to be evaluated, as in the seq function of the base package.}
\item{by}{Step size by which the defined parameter space is mapped out, as in the seq function of the base package.}
\item{data}{A vector of your data points.}
\item{priors}{A matrix of densities with in each column a density of a specific prior mapped on the paramater space that
is equal to the parameter space that is supplied using the from, to, by statements. E.g. the parameter space
runs from -10 to 10 in steps of 0.01 than your density of a standard normal distribution should be obtained
using dnorm(x = seq(from = -10, to = 10, by = 0.01), mean = 0, sd = 1). The first column will thus describe this
density using 2001 rows and all other columns should use the same density mapping to the parameter space.}
\item{lb.bench}{Lower boundary of the benchmark prior.}
\item{ub.bench}{Upper boundary of the benchmark prior.}
\item{n.iter}{The number of iterations that is used to obtain the posterior distribution of the data and the benchmark prior
note that only half of these iterations will be used to obtain samples, the other half is used for adaptation and
burnin.}
}

\references{
Bousquet, N. (2008) Diagnostics of prior-data agreement in applied bayesian analysis.
\emph{Journal of Applied Statistics}, \bold{35(9)}, 1011-1029.

Veen, D., Stoel, D., Schalken, N., van de Schoot, R. (2017). Using the Data Agreement Criterion to Rank Experts' Beliefs. \emph{arXiv preprint arXiv:1709.03736.}
}


\examples{
y <- Turnover$actual.data
from <- -10 
to <- 10
by <- .01
xx <- seq(from, to, by)
priors <- matrix(NA,ncol=4,nrow=length(xx))
for(i in 1:4){
  priors[,i] <- dnorm(xx, Turnover$expert.priors[i ,1], Turnover$expert.priors[i, 2])
}
out <- DAC.uniform(from = from, to = to, by = by, data = Turnover$actual.data, 
                   priors = priors ,lb.bench = 0, ub.bench = 5, n.iter = 4000)
out
}
