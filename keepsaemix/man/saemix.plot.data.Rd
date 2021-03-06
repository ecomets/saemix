\name{saemix.plot.data}
\alias{saemix.plot.data}
\alias{saemix.plot.convergence}
\alias{saemix.plot.llis}
\alias{saemix.plot.obsvspred}
\alias{saemix.plot.distribresiduals}
\alias{saemix.plot.scatterresiduals}
\alias{saemix.plot.fits}
\alias{saemix.plot.distpsi}
\alias{saemix.plot.randeff}
\alias{saemix.plot.correlations}
\alias{saemix.plot.parcov}
\alias{saemix.plot.randeffcov}
\alias{saemix.plot.npde}
\alias{saemix.plot.vpc}

\alias{saemix.plot.parcov.aux}
\alias{compute.sres}
\alias{compute.eta.map}

\title{ Functions implementing each type of plot in SAEM}
\description{ Several plots (selectable by the type argument) are currently available: convergence plot, individual plots, predictions versus observations, distribution plots, VPC, residual plots.
}
\usage{
saemix.plot.data(saemixObject, ...) 
saemix.plot.convergence(saemixObject,niter=0, ...) 
saemix.plot.llis(saemixObject, ...) 
saemix.plot.obsvspred(saemixObject, ...) 
saemix.plot.distribresiduals(saemixObject, ...) 
saemix.plot.scatterresiduals(saemixObject, ...) 
saemix.plot.fits(saemixObject, ...) 
saemix.plot.distpsi(saemixObject, ...) 
saemix.plot.randeff(saemixObject, ...) 
saemix.plot.correlations(saemixObject, ...) 
saemix.plot.parcov(saemixObject, ...) 
saemix.plot.randeffcov(saemixObject, ...) 
saemix.plot.npde(saemixObject, ...) 
saemix.plot.vpc(saemixObject,npc = FALSE, ...) 

saemix.plot.parcov.aux(saemixObject, partype = "p", ...)
compute.sres(saemixObject)
compute.eta.map(saemixObject)
}
\arguments{   
  \item{saemixObject}{ an object returned by the \code{\link{saemix}} function}
  \item{\dots}{ optional arguments passed to the plots}
  \item{npc}{ for VPC, computes Numerical Predictive Checks (currently not implemented)}
  \item{niter}{ the convergence plots are shown up to iteration "niter". Defaults to 0, which indicates the convergence plots should be plotted up to the maximal number of iterations set for the algorithm}
  \item{partype}{(this function is not user-level)}
}
\details{These functions implement plots different graphs related to the algorithm (convergence plots, likelihood estimation) as well as diagnostic graphs. A description is provided in the PDF documentation.

saemix.plot.parcov.aux, compute.sres and compute.eta.map are helper functions, not intended to be called by the user directly.

By default, the following plots are produced:
\describe{
  \item{saemix.plot.data:}{A spaghetti plot of the data, displaying the observed data y as a function of the regression variable (time for a PK application)}
  \item{saemix.plot.convergence:}{For each parameter in the model, this plot shows the evolution of the parameter estimate versus the iteration number}
  \item{saemix.plot.llis:}{Graph showing the evolution of the log-likelihood during the estimation by importance sampling}
  \item{saemix.plot.obsvspred:}{Plot of the predictions computed with the population parameters versus the observations (left), and plot of the predictions computed with the individual parameters versus the observations (right)}
  \item{saemix.plot.scatterresiduals:}{Scatterplot of the residuals versus the predictor (top) and versus predictions (bottom), for weighted residuals (population residuals, left), individual weighted residuals (middle) and npde (right).}
  \item{saemix.plot.distribresiduals:}{Distribution of the residuals, plotted as histogram (top) and as a QQ-plot (bottom), for weighted residuals (population residuals, left), individual weighted residuals (middle) and npde (right).}
  \item{saemix.plot.fits:}{Model fits. Individual fits are obtained using the individual parameters with the individual covariates. Population fits are obtained using the population parameters with the individual covariates (red) and the individual parameters with the individual covariates (green). By default the individual plots are displayed.}
  \item{saemix.plot.distpsi:}{Distribution of the parameters (conditional on covariates when some are included in the model). A histogram of individual parameter estimates can be overlayed on the plot, but it should be noted that the histogram does not make sense when there are covariates influencing the parameters and a warning will be displayed}
  \item{saemix.plot.randeff:}{Boxplot of the random effects}
  \item{saemix.plot.correlations:}{Correlation between the random effects}
  \item{saemix.plot.parcov:}{Plots of the estimates of the individual parameters versus the covariates, using scatterplot for continuous covariates, boxplot for categorical covariates}
  \item{saemix.plot.randeffcov:}{Plots of the estimates of the random effects versus the covariates, using scatterplot for continuous covariates, boxplot for categorical covariates}
  \item{saemix.plot.npde:}{Plots 4 graphs to evaluate the shape of the distribution of the normalised prediction distribution errors (npde)}
  \item{saemix.plot.vpc:}{Visual Predictive Check, with options to include the prediction intervals around the boundaries of the selected interval as well as around the median (50th percentile of the simulated data). Several methods are available to define binning on the X-axis (see methods in the PDF guide).}
}

Each plot can be customised by modifying options, either through a list of options set by the \code{\link{saemix.plot.setoptions}} function, or on the fly by passing an option in the call to the plot (see examples).
}
\value{None}
\references{ Kuhn E, Lavielle M. Maximum likelihood estimation in nonlinear mixed effects models. Computational Statistics and Data Analysis 49, 4 (2005), 1020-1038.

Comets E, Lavenu A, Lavielle M. SAEMIX, an R version of the SAEM algorithm. 20th meeting of the Population Approach Group in Europe, Athens, Greece (2011), Abstr 2173.
}
\author{ Emmanuelle Comets <emmanuelle.comets@inserm.fr>, Audrey Lavenu, Marc Lavielle.
}

\seealso{\code{\link{SaemixObject}},\code{\link{saemix}}, \code{\link{saemix.plot.setoptions}}, \code{\link{saemix.plot.select}}, \code{\link{plot.saemix}}}
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

saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE)

# Not run (strict time constraints for CRAN)
# saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)

# Simulate data and compute weighted residuals and npde
# saemix.fit<-compute.sres(saemix.fit)

# Data
# saemix.plot.data(saemix.fit)

# Convergence
# saemix.plot.convergence(saemix.fit)

# Individual plot for subject 1, smoothed
# saemix.plot.fits(saemix.fit,ilist=1,smooth=TRUE)

# Individual plot for subject 1 to 12, with ask set to TRUE 
# (the system will pause before a new graph is produced)
# saemix.plot.fits(saemix.fit,ilist=c(1:12),ask=TRUE)

# Diagnostic plot: observations versus population predictions
# par(mfrow=c(1,1))
# saemix.plot.obsvspred(saemix.fit,level=0,new=FALSE)

# LL by Importance Sampling
# saemix.plot.llis(saemix.fit)

# Scatter plot of residuals
# saemix.plot.scatterresiduals(saemix.fit)

# Boxplot of random effects
# saemix.plot.randeff(saemix.fit)

# Relationships between parameters and covariates
# saemix.plot.parcov(saemix.fit)

# Relationships between parameters and covariates, on the same page
# par(mfrow=c(3,2))
# saemix.plot.parcov(saemix.fit,new=FALSE)

# VPC, default options (10 bins, equal number of observations in each bin)
# Not run (time constraints for CRAN)
# saemix.plot.vpc(saemix.fit)

# VPC, user-defined breaks for binning
# Not run (time constraints for CRAN)
# saemix.plot.vpc(saemix.fit,vpc.method="user", vpc.breaks=c(0.4,0.8,1.5,2.5,4,5.5,8,10,13))
}
\keyword{plot}
