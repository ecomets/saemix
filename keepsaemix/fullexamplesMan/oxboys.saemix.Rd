\name{oxboys.saemix}
\docType{data}
\alias{oxboys.saemix}
\title{Heights of Boys in Oxford}
\description{The \code{oxboys.saemix} data frame has 234 rows and 4 columns.}
\usage{oxboys.saemix}
\format{       This data frame contains the following columns:
  \describe{     
  \item{Subject:}{  an ordered factor giving a unique identifier for each boy in the experiment
    }
    \item{age:}{  a numeric vector giving the standardized age (dimensionless)
    }
    \item{height:}{  a numeric vector giving the height of the boy (cm)
    }
    \item{Occasion:}{  an ordered factor - the result of converting 'age' from a continuous variable to a count so these slightly unbalanced data can be analyzed as balanced.
    }
  }
}
\details{      
These data are described in Goldstein (1987) as data on the height of a selection of boys from Oxford, England versus a standardized age. The dataset can be found in the package \code{nlme}.

We use an linear model for this data:
       y_ij = Base_i + slope_i x_ij +epsilon_ij 
}
\source{
Pinheiro, J. C. and Bates, D. M. (2000), _Mixed-Effects Models in S and S-PLUS_, Springer, New York.  (Appendix A.19)
}
\examples{ 

data(oxboys.saemix)
saemix.data<-saemixData(name.data=oxboys.saemix,header=TRUE,
  name.group=c("Subject"),name.predictors=c("age"),name.response=c("height"),
  units=list(x="yr",y="cm"))

growth.linear<-function(psi,id,xidep) {
# input:
#   psi : matrix of parameters (2 columns, base and slope)
#   id : vector of indices 
#   xidep : dependent variables (same nb of rows as length of id)
# returns:
#   a vector of predictions of length equal to length of id
  x<-xidep[,1]
  base<-psi[id,1]
  slope<-psi[id,2]
  f<-base+slope*x
  return(f)
}
saemix.model<-saemixModel(model=growth.linear,description="Linear model",
  psi0=matrix(c(140,1),ncol=2,byrow=TRUE,dimnames=list(NULL,c("base","slope"))),
  transform.par=c(1,0),covariance.model=matrix(c(1,1,1,1),ncol=2,byrow=TRUE), 
  error.model="constant")

saemix.options<-list(algorithms=c(1,1,1),nb.chains=1,seed=201004,
  save=FALSE,save.graphs=FALSE)

# plot the data
plot(saemix.data)

saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)

}
\keyword{datasets}
