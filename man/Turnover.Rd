\name{Turnover}
\alias{Turnover}
\docType{data}
\title{
Expert predictions and actual turnover for the first quarter of 2016
}
\description{
Data set containing predictions of four experts concerning the turnover for the first quarter of 2016 and the actual realisation of the turnover. The scale is transformed so as to protect business sensitive information..
}
\usage{data(Turnover)}
\format{
  A list of 2 providing the predictions of the experts and the actual data.
  \describe{
    \item{\code{expert.priors}}{Predictions of turnover by experts expressed in distributional form for the dsnorm function of the fGarch package. Each row representens one expert, column the mean, column two the sd and column three the skewness parameter.}
    \item{\code{actual.data}}{The realisation of turnover for the company.}
    }
}

\references{Veen, D., Stoel, D., Schalken, N., van de Schoot, R. (2017). \emph{Using the Data Agreement Criterion to Rank Experts' Beliefs}. Manuscript submitted for publication.
}

\examples{
head(Turnover)
}
