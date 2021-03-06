\name{map.saemix}
\alias{map.saemix}
\title{ Estimates of the individual parameters (conditional mode)}
\description{ Compute the estimates of the individual parameters PSI_i (conditional mode - Maximum A Posteriori)
}
\usage{ 
map.saemix(saemixObject)
}
\arguments{   
  \item{saemixObject}{ an object returned by the \code{\link{saemix}} function}
}
\details{ The MCMC procedure is used to estimate the conditional mode (or Maximum A Posteriori) m(phi_i |yi ; hat{theta}) = Argmax_{phi_i} p(phi_i |yi ; hat{theta})
 }
\value{ \item{saemixObject:}{returns the object with the estimates of the MAP parameters (see example for usage)}
}
\references{ Kuhn E, Lavielle M. Maximum likelihood estimation in nonlinear mixed effects models. Computational Statistics and Data Analysis 49, 4 (2005), 1020-1038.

Comets E, Lavenu A, Lavielle M. SAEMIX, an R version of the SAEM algorithm. 20th meeting of the Population Approach Group in Europe, Athens, Greece (2011), Abstr 2173.
}
\author{ Emmanuelle Comets <emmanuelle.comets@inserm.fr>, Audrey Lavenu, Marc Lavielle.
}
\seealso{\code{\link{SaemixObject}},\code{\link{saemix}}}
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

saemix.options<-list(algorithm=c(1,0,0),seed=632545,
  save=FALSE,save.graphs=FALSE)

saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)

# Estimating the individual parameters using the result of saemix 
# & returning the result in the same object

saemix.fit<-map.saemix(saemix.fit)
}
\keyword{models}
