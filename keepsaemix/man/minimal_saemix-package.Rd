\name{saemix-package}
\alias{saemix-package}

\docType{package}
\title{ Stochastic Approximation Expectation Maximization (SAEM) algorithm for non-linear mixed effects models}
\description{
\tabular{ll}{
  - Computing the maximum likelihood estimator of the population parameters, without any approximation \cr of the model (linearization, quadrature approximation, . . . ), using the Stochastic Approximation \cr Expectation Maximization (SAEM) algorithm \cr
  - Estimation of the Fisher Information matrix \cr
  - Estimation of the individual parameters \cr
  - Estimation of the likelihood \cr
  - Plot convergence graphs \cr
  }
}
\details{ \tabular{ll}{ Package: \tab saemix\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2015-09-01\cr
License: \tab GPL (>=) 1.2\cr
LazyLoad: \tab yes\cr
}

The SAEM package includes a number of undocumented functions, which are not meant to be used directly by the user.
\describe{
  \item{default}{setdefault}
  \item{computational functions}{cutoff,cutoff.max, cutoff.eps, cutoff.res, compute.Uy, compute.Uy.nocov, conditional.distribution, gqg.mlx}
  \item{distributions}{normcdf, norminv}
  \item{error model}{error}
  \item{sampling}{trnd.mlx, tpdf.mlx, gammarnd.mlx}
  \item{parameter transformations}{transpsi, transphi, dtransphi}
}
}
\author{ Emmanuelle Comets <emmanuelle.comets@inserm.fr>, Audrey Lavenu, Marc Lavielle.
}
\references{ Kuhn, E., and Lavielle, M. Maximum likelihood estimation in nonlinear mixed effects models. Computational Statistics and Data Analysis 49, 4 (2005), 1020-1038.

}
\keyword{saemix}
\keyword{model}
\seealso{nlme,\code{\link{SaemixData}},\code{\link{SaemixModel}}, \code{\link{SaemixObject}},\code{\link{saemix}}}
\examples{
data(theo.saemix)

saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, 
  name.group=c("Id"),name.predictors=c("Dose","Time"),
  name.response=c("Concentration"),name.covariates=c("Weight","Sex"),
  units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")

plot(saemix.data)
}
