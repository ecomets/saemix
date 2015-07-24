\name{conddist.saemix}
\alias{conddist.saemix}
\title{Estimate conditional mean and variance of individual parameters using the MCMC algorithm}
\description{When the parameters of the model have been estimated, we can estimate the individual parameters (psi_i). 

Let hat{theta} be the estimated value of theta computed with the SAEM algorithm and let p(phi_i |y_i; hat{theta}) be the conditional distribution of phi_i for 1<=i<=N . We use the MCMC procedure used in the SAEM algorithm to estimate these conditional distributions. We empirically estimate the conditional mean E(phi_i |y_i; hat{theta}) and the conditional standard deviation sd(phi_i |y_i; hat{theta}).
}
\usage{
conddist.saemix(saemixObject,nsamp=1,max.iter=NULL,...)
}
\arguments{
  \item{saemixObject}{ an object returned by the \code{\link{saemix}} function}
  \item{nsamp}{Number of samples to be drawn in the conditional distribution for each subject. Defaults to 1}
  \item{max.iter}{Maximum number of iterations for the computation of the conditional estimates. Defaults to twice the total number of iterations (sum(saemixObject["options"]$nbiter.saemix)*2)}
  \item{\dots}{ optional arguments passed to the plots. Plots will appear if the option displayProgress in the saemixObject object is TRUE}
}
\details{See PDF documentation for details of the computation. Briefly, the MCMC algorithm is used to obtain samples from the individual conditional distributions of the parameters. The algorithm is initialised for each subject to the conditional estimate of the individual parameters obtained at the end of the SAEMIX fit. A convergence criterion is used to ensure convergence of the mean and variance of the conditional distributions. When nsamp>1, several chains of the MCMC algorithm are run in parallel to obtain samples from the conditional distributions, and the convergence criterion must be achieved for all chains. When nsamp>1, the estimate of the conditional mean is obtained by averaging over the different samples.

The shrinkage for any given parameter for the conditional estimate is obtained as 
  
  Sh=1-var(eta_i)/omega(eta)

where var(eta_i) is the empirical variance of the estimates of the individual random effects, and omega(eta) is the estimated variance.

  The function adds or modifies the following elements in the results:
\describe{
  \item{cond.mean.phi}{Conditional mean of the individual distribution of the parameters (obtained as the mean of the samples)}
  \item{cond.var.phi}{Conditional variance of the individual distribution of the parameters (obtained as the mean of the estimated variance of the samples)}
  \item{cond.shrinkage}{Estimate of the shrinkage for the conditional estimates}
  \item{cond.mean.eta}{Conditional mean of the individual distribution of the parameters (obtained as the mean of the samples)}
  \item{phi.samp}{An array with 3 dimensions, giving nsamp samples from the conditional distributions of the individual parameters}
  \item{phi.samp.var}{The estimated individual variances for the sampled parameters phi.samp}
}
A warning is output if the maximum number of iterations is reached without convergence (the maximum number of iterations is saemix.options$nbiter.saemix[2]).
}
\references{ Kuhn E, Lavielle M. Maximum likelihood estimation in nonlinear mixed effects models. Computational Statistics and Data Analysis 49, 4 (2005), 1020-1038.

Comets E, Lavenu A, Lavielle M. SAEMIX, an R version of the SAEM algorithm. 20th meeting of the Population Approach Group in Europe, Athens, Greece (2011), Abstr 2173.
}
\author{ Emmanuelle Comets <emmanuelle.comets@inserm.fr>, Audrey Lavenu, Marc Lavielle.
}
\seealso{\code{\link{SaemixData}},\code{\link{SaemixModel}}, \code{\link{SaemixObject}},\code{\link{saemixControl}},\code{\link{saemix}}}
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
# saemix.fit<-conddist.saemix(saemix.fit,nsamp=3)

# First sample from the conditional distribution 
# (a N (nb of subject) by nb.etas (nb of parameters) matrix)
# saemix.fit["results"]["phi.samp"][,,1]

# Second sample
# saemix.fit["results"]["phi.samp"][,,2]

}
\keyword{model}
