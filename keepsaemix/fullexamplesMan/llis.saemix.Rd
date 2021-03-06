\name{llis.saemix}
\alias{llis.saemix}
\title{ Log-likelihood using Importance Sampling}
\description{ Estimate the log-likelihood using Importance Sampling
}
\usage{
llis.saemix(saemixObject)
}
\arguments{
  \item{saemixObject}{ an object returned by the \code{\link{saemix}} function}
}
\details{The likelihood of the observations is estimated without any approximation using a Monte-Carlo approach (see documentation).}
\value{ the log-likelihood estimated by Importance Sampling
}
\references{ Kuhn E, Lavielle M. Maximum likelihood estimation in nonlinear mixed effects models. Computational Statistics and Data Analysis 49, 4 (2005), 1020-1038.

Comets E, Lavenu A, Lavielle M. SAEMIX, an R version of the SAEM algorithm. 20th meeting of the Population Approach Group in Europe, Athens, Greece (2011), Abstr 2173.
}
\author{ Emmanuelle Comets <emmanuelle.comets@inserm.fr>, Audrey Lavenu, Marc Lavielle.
}

\seealso{\code{\link{SaemixObject}},\code{\link{saemix}},\code{\link{llgq.saemix}} }
\examples{ 
# Running the main algorithm to estimate the population parameters
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

# Not run (time constraints for CRAN)
# saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)

# Estimating the likelihood by importance sampling using the result of saemix 
# & returning the result in the same object
# Not run (time constraints for CRAN)
# saemix.fit<-llis.saemix(saemix.fit)
}
\keyword{models}
