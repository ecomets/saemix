\name{saemix-package}
\alias{saemix-package}
\alias{cutoff}
\alias{cutoff.max}
\alias{cutoff.eps}
\alias{cutoff.res}
\alias{error}
\alias{mydiag}
\alias{transpsi}
\alias{transphi}
\alias{dtransphi}
\alias{normcdf}
\alias{norminv}
\alias{compute.Uy}
\alias{conditional.distribution}
\alias{trnd.mlx}
\alias{tpdf.mlx}
\alias{gammarnd.mlx}
\alias{gqg.mlx}
\alias{derivphi}

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
Version: \tab 0.9\cr
Date: \tab 2010-09-19\cr
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
Monolix32_UsersGuide.pdf (http://software.monolix.org/sdoms/software/)
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

model1cpt<-function(psi,id,xidep) { 
	  dose<-xidep[,1]
	  tim<-xidep[,2]  
	  ka<-psi[id,1]
	  V<-psi[id,2]
	  CL<-psi[id,3]
	  k<-CL/V
	  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
	  return(ypred)
}

saemix.model<-saemixModel(model=model1cpt,
  description="One-compartment model with first-order absorption",  
  psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3, byrow=TRUE,dimnames=list(NULL, 
  c("ka","V","CL"))),transform.par=c(1,1,1), 
  covariate.model=matrix(c(0,1,0,0,0,0),ncol=3,byrow=TRUE),fixed.estim=c(1,1,1),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE), 
  omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE), error.model="constant")

saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE)

# Not run (strict time constraints for CRAN)
# saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)

# print(saemix.fit)
# plot(saemix.fit)
}
