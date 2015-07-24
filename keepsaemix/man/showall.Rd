\name{showall}
\alias{showall}
\title{ Prints out an extensive summary of an object}
\description{ This function shows an extensive summary of an object, and is used mainly to visualise the majority of the elements of an object}
\usage{
showall(object)
}
\arguments{
  \item{object}{ showall methods are available for objects of typeSaemixData, SaemixModel and SaemixObject}
}
\details{ None}
\value{None}
\seealso{\code{\link{SaemixData}},\code{\link{SaemixModel}}, \code{\link{SaemixObject}}}
\examples{
# A SaemixData object
data(theo.saemix)
saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, 
  name.group=c("Id"),name.predictors=c("Dose","Time"),
  name.response=c("Concentration"),name.covariates=c("Weight","Sex"),
  units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
showall(saemix.data)

# A SaemixModel object
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
showall(saemix.model)
}
\keyword{methods}
