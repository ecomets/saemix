

#' Methods for Function coef
#' 
#' Methods for function \code{coef}
#' 
#' 
#' @name coef-methods
#' @aliases coef-methods coef,ANY-method
#' @docType methods
#' @section Methods: \describe{
#' 
#' \item{list("signature(x = \"ANY\")")}{ default coef function ?}
#' 
#' \item{list("signature(x = \"SaemixObject\")")}{ extracts coefficients from
#' an SaemixObject}
#' 
#' }
#' @keywords methods
NULL





#' Methods for Function initialize
#' 
#' Constructor functions for Classes in the saemix package (not user-level)
#' 
#' 
#' @name initialize-methods
#' @docType methods
#' @section Methods: \describe{
#' 
#' \item{list("signature(.Object = \"SaemixData\")")}{ create a SaemixData
#' object. Please use the \code{\link{saemixData}} function.}
#' 
#' \item{list("signature(.Object = \"SaemixModel\")")}{ create a SaemixModel
#' object Please use the \code{\link{saemixModel}} function.}
#' 
#' \item{list("signature(.Object = \"SaemixObject\")")}{ create a SaemixObject
#' object. This object is obtained after a successful call to
#' \code{\link{saemix}}}
#' 
#' \item{list("signature(.Object = \"SaemixRepData\")")}{ create a
#' SaemixRepData object}
#' 
#' \item{list("signature(.Object = \"SaemixRes\")")}{ create a SaemixRes
#' object}
#' 
#' \item{list("signature(.Object = \"SaemixSimData\")")}{ create a
#' SaemixSimData object} }
#' @keywords methods
NULL





#' Methods for Function plot
#' 
#' Methods for function \code{plot}
#' 
#' 
#' @name plot-methods
#' @aliases plot-methods plot,ANY-method
#' @docType methods
#' @section Methods: \describe{
#' 
#' \item{list("signature(x = \"ANY\")")}{ default plot function ?}
#' 
#' \item{list("signature(x = \"SaemixData\")")}{ Plots the data. Defaults to a
#' spaghetti plot of response versus predictor, with lines joining the data for
#' one individual.}
#' 
#' \item{list("signature(x = \"SaemixModel\")")}{ Plots prediction of the model
#' }
#' 
#' \item{list("signature(x = \"SaemixObject\")")}{ This method gives access to
#' a number of plots that can be performed on a SaemixObject}
#' 
#' \item{list("signature(x = \"SaemixSimData\")")}{ Plots simulated datasets} }
#' @keywords methods
NULL





#' Methods for Function predict
#' 
#' Methods for function \code{predict}
#' 
#' 
#' @name predict-methods
#' @aliases predict-methods predict,ANY-method
#' @docType methods
#' @section Methods: \describe{
#' 
#' \item{list("signature(object = \"ANY\")")}{Default predict functions}
#' 
#' \item{list("signature(object = \"SaemixObject\")")}{Computes predictions
#' using the results of an SAEM fit} }
#' @keywords methods
NULL





#' Methods for Function print
#' 
#' Prints a summary of an object
#' 
#' 
#' @name print-methods
#' @aliases print.saemix print-methods print,ANY-method
#' @docType methods
#' @section Methods: \describe{
#' 
#' \item{list("signature(x = \"ANY\")")}{ Default print function }
#' 
#' \item{list("signature(x = \"SaemixData\")")}{ Prints a summary of a
#' SaemixData object }
#' 
#' \item{list("signature(x = \"SaemixModel\")")}{ Prints a summary of a
#' SaemixModel object }
#' 
#' \item{list("signature(x = \"SaemixObject\")")}{ Prints a summary of the
#' results from a SAEMIX fit }
#' 
#' \item{list("signature(x = \"SaemixRes\")")}{ Not user-level } }
#' @keywords methods
NULL





#' Methods for Functions psi, phi and eta
#' 
#' These methods are used to access estimates of individual parameters and
#' random effects
#' 
#' 
#' @name psi-methods
#' @aliases psi-methods phi-methods eta-methods psi,SaemixObject-method
#' phi,SaemixObject-method eta,SaemixObject-method
#' @docType methods
#' @section Methods: \describe{
#' 
#' \item{list("signature(object = \"SaemixObject\")")}{ please refer to the PDF
#' documentation for the models} }
#' @keywords methods
NULL





#' Class "SaemixData"
#' 
#' An object of the SaemixData class, representing a longitudinal data
#' structure, used by the SAEM algorithm.
#' 
#' 
#' @name SaemixData-class
#' @aliases SaemixData-class SaemixData [<-,SaemixData-method
#' [,SaemixData-method initialize,SaemixData-method plot,SaemixData-method
#' print,SaemixData-method showall,SaemixData-method show,SaemixData-method
#' plot,SaemixData print,SaemixData showall,SaemixData show,SaemixData
#' summary,SaemixData summary,SaemixData-method read.saemixData
#' read.saemixData,SaemixData-method SaemixRepData-class SaemixRepData
#' [<-,SaemixRepData-method [,SaemixRepData-method
#' initialize,SaemixRepData-method show,SaemixRepData-method
#' SaemixSimData-class SaemixSimData [<-,SaemixSimData-method
#' [,SaemixSimData-method initialize,SaemixSimData-method
#' plot,SaemixSimData-method show,SaemixSimData-method
#' @docType class
#' @section Objects from the Class: An object of the SaemixData class can be
#' created by using the function \code{\link{saemixData}}
#' @author Emmanuelle Comets <emmanuelle.comets@@inserm.fr>, Audrey Lavenu,
#' Marc Lavielle.
#' @seealso \code{\link{saemixData}},\code{\link{SaemixModel}},
#' \code{\link{saemixControl}},\code{\link{saemix}}
#' @references Kuhn E, Lavielle M. Maximum likelihood estimation in nonlinear
#' mixed effects models. Computational Statistics and Data Analysis 49, 4
#' (2005), 1020-1038.
#' 
#' Comets E, Lavenu A, Lavielle M. SAEMIX, an R version of the SAEM algorithm.
#' 20th meeting of the Population Approach Group in Europe, Athens, Greece
#' (2011), Abstr 2173.
#' @keywords classes
#' @examples
#' 
#' showClass("SaemixData")
#' 
#' # Specifying column names
#' data(theo.saemix)
#' saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, 
#'   name.group=c("Id"),name.predictors=c("Dose","Time"),
#'   name.response=c("Concentration"),name.covariates=c("Weight","Sex"),
#'   units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
#' 
#' # Specifying column numbers
#' data(theo.saemix)
#' saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA,
#'   name.group=1,name.predictors=c(2,3),name.response=c(4), name.covariates=5:6, 
#'   units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
#' 
#' # No column names specified, using automatic recognition of column names
#' data(PD1.saemix)
#' saemix.data<-saemixData(name.data=PD1.saemix,header=TRUE, 
#'   name.covariates=c("gender"),units=list(x="mg",y="-",covariates=c("-")))
#' 
NULL





#' Internal saemix objects
#' 
#' Internal saemix objects.
#' 
#' These are not to be called by the user.
#' 
#' @aliases .First.lib plotnpde
#' @keywords internal
NULL





#' Class "SaemixModel"
#' 
#' An object of the SaemixModel class, representing a non-linear mixed-effect
#' model structure, used by the SAEM algorithm.
#' 
#' 
#' @name SaemixModel-class
#' @aliases SaemixModel-class SaemixModel [<-,SaemixModel-method
#' [,SaemixModel-method initialize,SaemixModel-method plot,SaemixModel-method
#' print,SaemixModel-method showall,SaemixModel-method show,SaemixModel-method
#' plot,SaemixModel print,SaemixModel showall,SaemixModel show,SaemixModel
#' summary,SaemixModel summary,SaemixModel-method
#' @docType class
#' @section Objects from the Class: An object of the SaemixModel class can be
#' created by using the function \code{\link{saemixModel}}
#' @author Emmanuelle Comets <emmanuelle.comets@@inserm.fr>, Audrey Lavenu,
#' Marc Lavielle.
#' @seealso \code{\link{SaemixData}},\code{\link{saemixModel}},
#' \code{\link{saemixControl}},\code{\link{saemix}}
#' @references Kuhn E, Lavielle M. Maximum likelihood estimation in nonlinear
#' mixed effects models. Computational Statistics and Data Analysis 49, 4
#' (2005), 1020-1038.
#' 
#' Comets E, Lavenu A, Lavielle M. SAEMIX, an R version of the SAEM algorithm.
#' 20th meeting of the Population Approach Group in Europe, Athens, Greece
#' (2011), Abstr 2173.
#' @keywords classes
#' @examples
#' 
#' showClass("SaemixModel")
#' 
NULL





#' Class "SaemixObject"
#' 
#' An object of the SaemixObject class, storing the results obtained by a call
#' to the SAEM algorithm
#' 
#' Details of the algorithm can be found in the pdf file accompanying the
#' package.
#' 
#' @name SaemixObject-class
#' @aliases SaemixObject-class SaemixObject [<-,SaemixObject-method
#' [,SaemixObject-method initialize,SaemixObject-method
#' plot,SaemixObject-method print,SaemixObject-method
#' predict,SaemixObject-method showall,SaemixObject-method
#' show,SaemixObject-method print,SaemixObject predict,SaemixObject
#' showall,SaemixObject show,SaemixObject summary,SaemixObject
#' summary,SaemixObject-method coef,SaemixObject coef,SaemixObject-method
#' SaemixRes-class SaemixRes [<-,SaemixRes-method [,SaemixRes-method
#' initialize,SaemixRes-method print,SaemixRes-method showall,SaemixRes-method
#' show,SaemixRes-method print,SaemixRes showall,SaemixRes show,SaemixRes
#' @docType class
#' @section Objects from the Class: Objects are created by a call
#' \code{saemix()}.
#' @author Emmanuelle Comets <emmanuelle.comets@@inserm.fr>, Audrey Lavenu,
#' Marc Lavielle.
#' @seealso \code{\link{SaemixData}},\code{\link{SaemixModel}},
#' \code{\link{saemixControl}},\code{\link{saemix}},\code{\link{plot.saemix}},
#' \code{\link{saemix.plot.data}}
#' @references Kuhn E, Lavielle M. Maximum likelihood estimation in nonlinear
#' mixed effects models. Computational Statistics and Data Analysis 49, 4
#' (2005), 1020-1038.
#' 
#' Comets E, Lavenu A, Lavielle M. SAEMIX, an R version of the SAEM algorithm.
#' 20th meeting of the Population Approach Group in Europe, Athens, Greece
#' (2011), Abstr 2173.
#' @keywords classes
#' @examples
#' 
#' showClass("SaemixObject")
#' 
NULL





#' General plot function from SAEM
#' 
#' Several plots (selectable by the type argument) are currently available:
#' convergence plot, individual plots, predictions versus observations,
#' distribution plots, VPC, residual plots.
#' 
#' This is the generic plot function for an SaemixObject object, which
#' implements different graphs related to the algorithm (convergence plots,
#' likelihood estimation) as well as diagnostic graphs. A description is
#' provided in the PDF documentation. Arguments such as main, xlab, etc... that
#' can be given to the generic plot function may be used, and will be
#' interpreted according to the type of plot that is to be drawn.
#' 
#' A special argument plot.type can be set to determine the type of plot; it
#' can be one of: \describe{ \item{data:}{A spaghetti plot of the data,
#' displaying the observed data y as a function of the regression variable
#' (time for a PK application)} \item{convergence:}{For each parameter in the
#' model, this plot shows the evolution of the parameter estimate versus the
#' iteration number} \item{likelihood:}{Graph showing the evolution of the
#' log-likelihood during the estimation by importance sampling}
#' \item{observations.vs.predictions:}{Plot of the predictions computed with
#' the population parameters versus the observations (left), and plot of the
#' predictions computed with the individual parameters versus the observations
#' (right)} \item{residuals.scatter:}{Scatterplot of the residuals versus the
#' predictor (top) and versus predictions (bottom), for weighted residuals
#' (population residuals, left), individual weighted residuals (middle) and
#' npde (right).} \item{residuals.distribution:}{Distribution of the residuals,
#' plotted as histogram (top) and as a QQ-plot (bottom), for weighted residuals
#' (population residuals, left), individual weighted residuals (middle) and
#' npde (right).} \item{individual.fit:}{Individual fits are obtained using the
#' individual parameters with the individual covariates}
#' \item{population.fit:}{Population fits are obtained using the population
#' parameters with the individual covariates} \item{both.fit:}{Individual fits,
#' superposing fits obtained using the population parameters with the
#' individual covariates (red) and using the individual parameters with the
#' individual covariates (green)} \item{marginal.distribution:}{Distribution of
#' the parameters (conditional on covariates when some are included in the
#' model). A histogram of individual parameter estimates can be overlayed on
#' the plot, but it should be noted that the histogram does not make sense when
#' there are covariates influencing the parameters and a warning will be
#' displayed} \item{random.effects:}{Boxplot of the random effects}
#' \item{correlations:}{Correlation between the random effects}
#' \item{parameters.vs.covariates:}{Plots of the estimates of the individual
#' parameters versus the covariates, using scatterplot for continuous
#' covariates, boxplot for categorical covariates}
#' \item{randeff.vs.covariates:}{Plots of the estimates of the random effects
#' versus the covariates, using scatterplot for continuous covariates, boxplot
#' for categorical covariates} \item{npde:}{Plots 4 graphs to evaluate the
#' shape of the distribution of the normalised prediction distribution errors
#' (npde)} \item{vpc:}{Visual Predictive Check, with options to include the
#' prediction intervals around the boundaries of the selected interval as well
#' as around the median (50th percentile of the simulated data).} } In
#' addition, the following values for plot.type produce a series of plots:
#' \describe{ \item{reduced:}{ produces the following plots: plot of the data,
#' convergence plots, plot of the likelihood by importance sampling (if
#' computed), plots of observations versus predictions. This is the default
#' behaviour of the plot function applied to an SaemixObject object}
#' \item{full:}{ produces the following plots: plot of the data, convergence
#' plots, plot of the likelihood by importance sampling (if computed), plots of
#' observations versus predictions, scatterplots and distribution of residuals,
#' VPC, npde, boxplot of the random effects, distribution of the parameters,
#' correlations between random effects, plots of the relationships between
#' individually estimated parameters and covariates, plots of the relationships
#' between individually estimated random effects and covariates}
#' 
#' Each plot can be customised by modifying options, either through a list of
#' options set by the \code{\link{saemix.plot.setoptions}} function, or on the
#' fly by passing an option in the call to the plot (see examples). }
#' 
#' @aliases plot.saemix plot,SaemixObject plot
#' @param x an object returned by the \code{\link{saemix}} function
#' @param y empty
#' @param \dots optional arguments passed to the plots
#' @return None
#' @author Emmanuelle Comets <emmanuelle.comets@@inserm.fr>, Audrey Lavenu,
#' Marc Lavielle.
#' @seealso \code{\link{SaemixObject}},\code{\link{saemix}},
#' \code{\link{saemix.plot.setoptions}}, \code{\link{saemix.plot.select}},
#' \code{\link{saemix.plot.data}}
#' @references Kuhn E, Lavielle M. Maximum likelihood estimation in nonlinear
#' mixed effects models. Computational Statistics and Data Analysis 49, 4
#' (2005), 1020-1038.
#' 
#' Comets E, Lavenu A, Lavielle M. SAEMIX, an R version of the SAEM algorithm.
#' 20th meeting of the Population Approach Group in Europe, Athens, Greece
#' (2011), Abstr 2173.
#' @keywords plot
#' @examples
#' 
#' data(theo.saemix)
#' 
#' saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, 
#'   name.group=c("Id"),name.predictors=c("Dose","Time"),
#'   name.response=c("Concentration"),name.covariates=c("Weight","Sex"),
#'   units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
#' 
#' model1cpt<-function(psi,id,xidep) { 
#' 	  dose<-xidep[,1]
#' 	  tim<-xidep[,2]  
#' 	  ka<-psi[id,1]
#' 	  V<-psi[id,2]
#' 	  CL<-psi[id,3]
#' 	  k<-CL/V
#' 	  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
#' 	  return(ypred)
#' }
#' 
#' saemix.model<-saemixModel(model=model1cpt,
#'   description="One-compartment model with first-order absorption", 
#'   psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3, byrow=TRUE,
#'   dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1),
#'   covariate.model=matrix(c(0,1,0,0,0,0),ncol=3,byrow=TRUE),fixed.estim=c(1,1,1),
#'   covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
#'   omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="constant")
#' 
#' saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE)
#' 
#' # Not run (strict time constraints for CRAN)
#' # saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)
#' 
#' # Set of default plots
#' # plot(saemix.fit)
#' 
#' # Data
#' # plot(saemix.fit,plot.type="data")
#' 
#' # Convergence
#' # plot(saemix.fit,plot.type="convergence")
#' 
#' # Individual plot for subject 1, smoothed
#' # plot(saemix.fit,plot.type="individual.fit",ilist=1,smooth=TRUE)
#' 
#' # Individual plot for subject 1 to 12, with ask set to TRUE 
#' # (the system will pause before a new graph is produced)
#' # plot(saemix.fit,plot.type="individual.fit",ilist=c(1:12),ask=TRUE)
#' 
#' # Diagnostic plot: observations versus population predictions
#' # par(mfrow=c(1,1))
#' # plot(saemix.fit,plot.type="observations.vs.predictions",level=0,new=FALSE)
#' 
#' # LL by Importance Sampling
#' # plot(saemix.fit,plot.type="likelihood")
#' 
#' # Scatter plot of residuals
#' # Data will be simulated to compute weighted residuals and npde
#' # the results shall be silently added to the object saemix.fit
#' # plot(saemix.fit,plot.type="residuals.scatter")
#' 
#' # Boxplot of random effects
#' # plot(saemix.fit,plot.type="random.effects")
#' 
#' # Relationships between parameters and covariates
#' # plot(saemix.fit,plot.type="parameters.vs.covariates")
#' 
#' # Relationships between parameters and covariates, on the same page
#' # par(mfrow=c(3,2))
#' # plot(saemix.fit,plot.type="parameters.vs.covariates",new=FALSE)
#' 
#' # VPC
#' # Not run (time constraints for CRAN)
#' # plot(saemix.fit,plot.type="vpc")
#' 
NULL





#' Methods for Function showall
#' 
#' Prints an extensive summary of an object
#' 
#' 
#' @name showall-methods
#' @docType methods
#' @section Methods: \describe{ \item{list("signature(x = \"SaemixData\")")}{
#' Prints a extensive summary of a SaemixData object }
#' 
#' \item{list("signature(x = \"SaemixModel\")")}{ Prints a extensive summary of
#' a SaemixModel object }
#' 
#' \item{list("signature(x = \"SaemixObject\")")}{ Prints a extensive summary
#' of the results from a SAEMIX fit }
#' 
#' \item{list("signature(x = \"SaemixRes\")")}{ Not user-level } }
#' @keywords methods
NULL





#' Prints out an extensive summary of an object
#' 
#' This function shows an extensive summary of an object, and is used mainly to
#' visualise the majority of the elements of an object
#' 
#' None
#' 
#' @param object showall methods are available for objects of typeSaemixData,
#' SaemixModel and SaemixObject
#' @return None
#' @seealso \code{\link{SaemixData}},\code{\link{SaemixModel}},
#' \code{\link{SaemixObject}}
#' @keywords methods
#' @examples
#' 
#' # A SaemixData object
#' data(theo.saemix)
#' saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, 
#'   name.group=c("Id"),name.predictors=c("Dose","Time"),
#'   name.response=c("Concentration"),name.covariates=c("Weight","Sex"),
#'   units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")
#' showall(saemix.data)
#' 
#' # A SaemixModel object
#' model1cpt<-function(psi,id,xidep) { 
#' 	  dose<-xidep[,1]
#' 	  tim<-xidep[,2]  
#' 	  ka<-psi[id,1]
#' 	  V<-psi[id,2]
#' 	  CL<-psi[id,3]
#' 	  k<-CL/V
#' 	  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
#' 	  return(ypred)
#' }
#' saemix.model<-saemixModel(model=model1cpt,
#'   description="One-compartment model with first-order absorption", 
#'   psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3, byrow=TRUE,
#'   dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1),
#'   covariate.model=matrix(c(0,1,0,0,0,0),ncol=3,byrow=TRUE),fixed.estim=c(1,1,1),
#'   covariance.model=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),
#'   omega.init=matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,byrow=TRUE),error.model="constant")
#' showall(saemix.model)
#' 
NULL





#' Methods for Function show
#' 
#' Prints a short summary of an object
#' 
#' 
#' @name show-methods
#' @docType methods
#' @section Methods: \describe{
#' 
#' \item{list("signature(x = \"ANY\")")}{ Default show function }
#' 
#' \item{list("signature(x = \"SaemixData\")")}{ Prints a short summary of a
#' SaemixData object }
#' 
#' \item{list("signature(x = \"SaemixModel\")")}{ Prints a short summary of a
#' SaemixModel object }
#' 
#' \item{list("signature(x = \"SaemixObject\")")}{ Prints a short summary of
#' the results from a SAEMIX fit }
#' 
#' \item{list("signature(x = \"SaemixRes\")")}{ Not user-level }
#' 
#' \item{list("signature(object = \"SaemixRepData\")")}{ Prints a short summary
#' of a SaemixRepData object }
#' 
#' \item{list("signature(object = \"SaemixSimData\")")}{ Prints a short summary
#' of a SaemixSimData object } }
#' @keywords methods
NULL





#' Methods for subset
#' 
#' This method is used to select a subset of an SaemixData object
#' 
#' 
#' @name subset-methods
#' @aliases subset-methods subset,ANY-method subset.SaemixData
#' subset,SaemixData-method
#' @docType methods
#' @section Methods: \describe{ \item{list("signature(x = \"ANY\")")}{ Default
#' subset function }
#' 
#' \item{list("signature(x = \"SaemixData\")")}{ A subset of the dataset can be
#' selected (see \code{subset} )} }
#' @keywords methods
NULL





#' Methods for Function summary
#' 
#' Methods for function \code{summary}
#' 
#' 
#' @name summary-methods
#' @aliases summary-methods summary,ANY-method
#' @docType methods
#' @section Methods: \describe{
#' 
#' \item{list("signature(x = \"ANY\")")}{ default summary function ?}
#' 
#' \item{list("signature(x = \"SaemixData\")")}{ summary of the data }
#' 
#' \item{list("signature(x = \"SaemixModel\")")}{ summary of the model }
#' 
#' \item{list("signature(x = \"SaemixObject\")")}{ summary of an SaemixObject}
#' 
#' }
#' @keywords methods
NULL



