\name{cow.saemix}
\docType{data}
\alias{cow.saemix}
\title{Evolution of the weight of 560 cows, in SAEM format}
\description{ \code{cow.saemix} contains data from winter wheat experiments.}
\usage{cow.saemix}
\format{   This data frame contains the following columns:
  \describe{     
  \item{cow:}{  the id.
    }
    \item{time:}{  time (days).
    }
    \item{weight:}{  weight of the cow (kg).
    }
    \item{birthyear:}{  year of birth (between 1988 and 1998).
    }
    \item{twin:}{  existence of a twin (no=1, yes=2).
    }
    \item{birthrank:}{ the rank of birth (beetween 3 and 7).
    }
  }
}
\details{      
The data used in this example is the evolution of the weight (in kg) of 560 cows. We have 9 or 10 measurements per subject.

We use an exponential growth model for this data:
       y_ij = A_i (1- B_i exp( - K_i t_ij)) +epsilon_ij 
}
\examples{ 

data(cow.saemix)
saemix.data<-saemixData(name.data=cow.saemix,header=TRUE,name.group=c("cow"), 
  name.predictors=c("time"),name.response=c("weight"), 
  name.covariates=c("birthyear","twin","birthrank"), 
  units=list(x="days",y="kg",covariates=c("yr","-","-")))

growthcow<-function(psi,id,xidep) {
# input:
#   psi : matrix of parameters (3 columns, a, b, k)
#   id : vector of indices 
#   xidep : dependent variables (same nb of rows as length of id)
# returns:
#   a vector of predictions of length equal to length of id
  x<-xidep[,1]
  a<-psi[id,1]
  b<-psi[id,2]
  k<-psi[id,3]
  f<-a*(1-b*exp(-k*x))
  return(f)
}
saemix.model<-saemixModel(model=growthcow,
  description="Exponential growth model", 
  psi0=matrix(c(700,0.9,0.02,0,0,0),ncol=3,byrow=TRUE, 
  dimnames=list(NULL,c("A","B","k"))),transform.par=c(1,1,1),fixed.estim=c(1,1,1), 
  covariate.model=matrix(c(0,0,0),ncol=3,byrow=TRUE), 
  covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE), 
  omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="constant")

saemix.options<-list(algorithms=c(1,1,1),nb.chains=1,nbiter.saemix=c(200,100), 
  seed=4526,save=FALSE,save.graphs=FALSE)

# Plotting the data
plot(saemix.data,xlab="Time (day)",ylab="Weight of the cow (kg)")

# Not run (strict time constraints for CRAN)
# saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)

}
\keyword{datasets}
