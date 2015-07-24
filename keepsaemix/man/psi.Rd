\name{psi}
\alias{psi}
\alias{phi}
\alias{eta}
\title{ Functions to extract the individual estimates of the parameters and random effects}
\description{ These three functions are used to access the estimates of individual parameters and random effects.}
\usage{
psi(object, indiv.par)
phi(object, indiv.par)
eta(object, indiv.par)
}
\arguments{
  \item{object}{ an object returned by the \code{\link{saemix}} function}
  \item{indiv.par}{ a string giving the type of estimate to be used (one of "map", for the Maximum A Posteriori estimate, or "eap", for conditional estimate). Defaults to "map"}
}
\details{ The psi_i represent the individual parameter estimates. In the SAEM algorithm, these parameters are assumed to be a transformation of a Gaussian random vector phi_i, where the phi_i can be written as a function of the individual random effects (eta_i), the covariate matrix (C_i) and the vector of fixed effects (mu):

phi_i = C_i mu + eta_i

More details can be found in the PDF documentation.}
\value{ These functions are used to access and output the estimates of parameters and random effects. When the object passed to the function does not contain these estimates, they are automatically computed. The object is then returned (invisibly) with these estimates added to the results. }
\references{ Kuhn E, Lavielle M. Maximum likelihood estimation in nonlinear mixed effects models. Computational Statistics and Data Analysis 49, 4 (2005), 1020-1038.

Comets E, Lavenu A, Lavielle M. SAEMIX, an R version of the SAEM algorithm. 20th meeting of the Population Approach Group in Europe, Athens, Greece (2011), Abstr 2173.
}
\author{ Emmanuelle Comets <emmanuelle.comets@inserm.fr>, Audrey Lavenu, Marc Lavielle.
}
\seealso{\code{\link{SaemixData}},\code{\link{SaemixModel}}, \code{\link{SaemixObject}}, \code{\link{saemixControl}}, \code{\link{plot.saemix}}}
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
  psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3, byrow=TRUE,
  dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1),
  covariate.model=matrix(c(0,1,0,0,0,0),ncol=3,byrow=TRUE),fixed.estim=c(1,1,1),
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
  omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="constant")

saemix.options<-list(algorithm=c(1,0,0),seed=632545,save=FALSE,save.graphs=FALSE)

# Not run (strict time constraints for CRAN)
# saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)

# psi(saemix.fit)
# phi(saemix.fit)
# eta(saemix.fit,indiv.par="eap")
}
\keyword{models}
