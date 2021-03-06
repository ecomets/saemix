\name{theo.saemix}
\docType{data}
\alias{theo.saemix}
\title{Pharmacokinetics of theophylline, in SAEM format}
\description{   The \code{theo.saemix} data frame has 132 rows and 6 columns of data from
  an experiment on the pharmacokinetics of theophylline. A column with gender was
  added to the original data for demo purposes, and contains simulated data.}
\usage{theo.saemix}
\format{   This data frame contains the following columns:
  \describe{     
  \item{Id:}{       an ordered factor with levels \code{1}, \dots, \code{12}
      identifying the subject on whom the observation was made.  The
      ordering is by Time at which the observation was made.
    }
    \item{Dose:}{       dose of theophylline administered orally to the subject (mg/kg).
    }
    \item{Time:}{       time since drug administration when the sample was drawn (hr).
    }
    \item{Concentration:}{       theophylline concentration in the sample (mg/L).
    }
    \item{Weight:}{       weight of the subject (kg).
    }
    \item{Sex:}{       gender of the subject (0=men, 1=women).
    }
  }
}
\details{   Boeckmann, Sheiner and Beal (1994) report data from a study by
  Dr. Robert Upton of the kinetics of the anti-asthmatic drug
  theophylline.  Twelve subjects were given oral doses of theophylline
  then serum concentrations were measured at 11 time points over the
  next 25 hours. In the present package \emph{npde}, we removed the data
  at time 0.

  These data are analyzed in Davidian and Giltinan (1995) and Pinheiro
  and Bates (2000) using a two-compartment open pharmacokinetic model.
  
  These data are also available in the library \code{datasets} under the name
  \code{Theoph} in a slightly modified format and including the data at time 0.

}
\source{   Boeckmann, A. J., Sheiner, L. B. and Beal, S. L. (1994), \emph{NONMEM
    Users Guide: Part V}, NONMEM Project Group, University of
  California, San Francisco.
    
  Davidian, M. and Giltinan, D. M. (1995) \emph{Nonlinear Models for
    Repeated Measurement Data}, Chapman & Hall (section 5.5, p. 145 and
  section 6.6, p. 176)

  Pinheiro, J. C. and Bates, D. M. (2000) \emph{Mixed-effects Models in
    S and S-PLUS}, Springer (Appendix A.29)
}
\examples{ 
data(theo.saemix)
saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, 
  name.group=c("Id"),name.predictors=c("Dose","Time"),
  name.response=c("Concentration"),name.covariates=c("Weight","Sex"),
  units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")

plot(saemix.data)
}
\keyword{datasets}
