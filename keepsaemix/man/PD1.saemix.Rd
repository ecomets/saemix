\name{PD1.saemix}
\docType{data}
\alias{PD1.saemix}
\alias{PD2.saemix}
\title{Data simulated according to an Emax response model, in SAEM format}
\description{ \code{PD1.saemix} contains data from winter wheat experiments.}
\usage{
PD1.saemix

PD2.saemix}
\format{   This data frame contains the following columns:
  \describe{     
  \item{subject:}{  the site number
    }
    \item{dose:}{  simulated dose.
    }
    \item{response:}{  simulated response.
    }
    \item{gender:}{   gender (0 for male, 1 for female).
    }
  }
}
\details{  
These examples were used by P. Girard and F. Mentre for the symposium dedicated to Comparison of Algorithms Using Simulated Data Sets and Blind Analysis, that took place in Lyon, France, September 2004.

The dataset contains 100 individuals, each receiving 3 different doses:(0, 10, 90), (5, 25, 65) or (0, 20, 30). It was assumed that doses were given in a cross-over study with sufficient wash out period to avoid carry over. Responses (y_ij) were simulated with the following pharmacodynamic model:

       y_ij = E0_i + D_ij Emax_i/(D_ij + ED50_i) +epsilon_ij 

The individual parameters were simulated according to

       log (E0_i) = log (E0) + eta_i1
       log (Emax_i) = log (Emax) + eta_i2
       log (E50_i) = log (E50) +  beta w_i + eta_i3

PD1.saemix contains the data simulated with a gender effect, beta=0.3.

PD2.saemix contains the data simulated without a gender effect, beta=0.

}
\source{  Girard P., Mentre F. Comparison of Algorithms Using Simulated Data Sets and Blind Analysis workshop, Lyon, France, September 2004. 
}
\examples{ 

data(PD1.saemix)
saemix.data<-saemixData(name.data=PD1.saemix,header=TRUE,name.group=c("subject"),
  name.predictors=c("dose"),name.response=c("response"),
  name.covariates=c("gender"), units=list(x="mg",y="-",covariates=c("-")))

modelemax<-function(psi,id,xidep) {
# input:
#   psi : matrix of parameters (3 columns, E0, Emax, EC50)
#   id : vector of indices 
#   xidep : dependent variables (same nb of rows as length of id)
# returns:
#   a vector of predictions of length equal to length of id
  dose<-xidep[,1]
  e0<-psi[id,1]
  emax<-psi[id,2]
  e50<-psi[id,3]
  f<-e0+emax*dose/(e50+dose)
  return(f)
}

# Plotting the data
plot(saemix.data,main="Simulated data PD1")

# Not run (strict time constraints for CRAN)

# Compare models with and without covariates with LL by Importance Sampling
# SE not computed
model1<-saemixModel(model=modelemax,description="Emax growth model", 
  psi0=matrix(c(20,300,20,0,0,0),ncol=3,byrow=TRUE,dimnames=list(NULL,
  c("E0","Emax","EC50"))), transform.par=c(1,1,1),
  covariate.model=matrix(c(0,0,0), ncol=3,byrow=TRUE),fixed.estim=c(1,1,1))

model2<-saemixModel(model=modelemax,description="Emax growth model", 
  psi0=matrix(c(20,300,20,0,0,0),ncol=3,byrow=TRUE,dimnames=list(NULL, 
  c("E0","Emax","EC50"))), transform.par=c(1,1,1),
  covariate.model=matrix(c(0,0,1), ncol=3,byrow=TRUE),fixed.estim=c(1,1,1))

saemix.options<-list(algorithms=c(0,1,1),nb.chains=3,seed=765754, 
  nbiter.saemix=c(500,300),save=FALSE,save.graphs=FALSE)

# Not run (strict time constraints for CRAN)
# fit1<-saemix(model1,saemix.data,saemix.options)
# fit2<-saemix(model2,saemix.data,saemix.options)
# wstat<-(-2)*(fit1["results"]["ll.is"]-fit2["results"]["ll.is"])

# cat("LRT test for covariate effect on EC50: p-value=",1-pchisq(wstat,1),"\n")
}
\keyword{datasets}
