\name{compute.sres}
\alias{compute.sres}
\alias{compute.eta.map}
\title{Compute residuals for a SAEMIX fit
}
\description{
This function computes standardised residuals (WRES), prediction discrepancies (pd) and normalised prediction distribution errors (npde).
}
\usage{
compute.sres(saemixObject)
compute.eta.map(saemixObject)
}
\arguments{
  \item{saemixObject}{ an object returned by the \code{\link{saemix}} function}
}
\details{ Computation of standardised residuals 

Function compute.eta is a helper function used to compute individual random effects given individual parameters.
}
\value{ Adds the following elements to the saemixObject.
\item{pd}{Prediction discrepancies}
\item{npd}{Normalised prediction distribution errors}
\item{wres}{Population weighted residuals}
\item{sim.data}{Simulated data}
\item{ypred}{The population mean predictions of the model, obtained as E(y_i(sim_k)) over the k simulations ( E_hat{theta} (f(x_i, psi_i)) )}
}

\references{K. Brendel, E. Comets, C. Laffont, C. Laveille, and F. Mentr\'e.
Metrics for external model evaluation with an application to the population
pharmacokinetics of gliclazide. \emph{Pharmaceutical Research}, 23:2036--49,
2006.
Kuhn, E., and Lavielle, M. Maximum likelihood estimation in nonlinear mixed effects models. Computational Statistics and Data Analysis 49, 4 (2005), 1020-1038.

Monolix32_UsersGuide.pdf (\url{http://software.monolix.org/sdoms/software/})
}
\author{ Emmanuelle Comets <emmanuelle.comets@inserm.fr>, Audrey Lavenu, Marc Lavielle.
}

\seealso{\code{\link{saemixData}},\code{\link{saemixModel}}, \code{\link{saemixControl}},\code{\link{saemix}},\code{\link{predict.saemix}}, \code{\link{simul.saemix}}}
\examples{
  data(theo.saemix)
  saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, name.group=c("Id"),name.predictors=c("Dose","Time"),name.response=c("Concentration"), name.covariates=c("Weight","Sex"),units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")

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
  saemix.model<-saemixModel(model=model1cpt,description="One-compartment model with first-order absorption", psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3, byrow=TRUE,dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1), covariate.model=matrix(c(0,1,0,0,0,0),ncol=3,byrow=TRUE),fixed.estim=c(1,1,1),
covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE), omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE), error.model="constant")

  saemix.options<-list(seed=632545)

  saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)

  saemix.fit<-compute.sres(saemix.fit)

# Scatterplot of residuals
  plot(saemix.fit,type="residuals.scatter")

# VPC
  plot(saemix.fit,type="vpc") 
}
\keyword{models}
