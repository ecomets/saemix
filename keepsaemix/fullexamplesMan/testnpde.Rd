\name{testnpde}
\alias{testnpde}
\alias{kurtosis}
\alias{skewness}
\title{Tests for normalised prediction distribution errors}
\description{
  Performs tests for the normalised prediction distribution errors returned by
\code{npde}
}
\usage{
testnpde(npde)
}
\arguments{
  \item{npde}{the vector of prediction distribution errors}
}
\details{
  Given a vector of normalised prediction distribution errors (npde), this function compares the npde to the standardised normal distribution N(0,1) using a Wilcoxon test of the mean, a Fisher test of the variance, and a Shapiro-Wilks test for normality. A global test is also reported.
  
  The helper functions \code{kurtosis} and \code{skewness} are called to compute the kurtosis and skewness of the distribution of the npde.
}
\value{
  a list containing 4 components:
  \item{Wilcoxon test of mean=0}{compares the mean of the npde
  to 0 using a Wilcoxon test}
  \item{variance test       }{compares the variance of the npde
  to 1 using a Fisher test}
  \item{SW test of normality}{compares the npde to the normal
  distribution using a Shapiro-Wilks test}
  \item{global test         }{an adjusted p-value corresponding to the minimum
  of the 3 previous p-values multiplied by the number of tests (3), or 1 if
  this p-value is larger than 1.}
}
\references{K. Brendel, E. Comets, C. Laffont, C. Laveille, and F. Mentr\'e.
Metrics for external model evaluation with an application to the population
pharmacokinetics of gliclazide. \emph{Pharmaceutical Research}, 23:2036--49,
2006.}
\author{Emmanuelle Comets <emmanuelle.comets@inserm.fr>}

\seealso{ \code{\link{saemix}}, \code{\link{saemix.plot.npde}}}
\keyword{models}
