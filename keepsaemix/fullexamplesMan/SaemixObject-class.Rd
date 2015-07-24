\name{SaemixObject-class}
\Rdversion{1.1}
\docType{class}
\alias{SaemixObject-class}
\alias{SaemixObject}
\alias{[<-,SaemixObject-method}
\alias{[,SaemixObject-method}
\alias{initialize,SaemixObject-method}
\alias{plot,SaemixObject-method}
\alias{print,SaemixObject-method}
\alias{predict,SaemixObject-method}
\alias{showall,SaemixObject-method}
\alias{show,SaemixObject-method}
\alias{print,SaemixObject}
\alias{predict,SaemixObject}
\alias{showall,SaemixObject}
\alias{show,SaemixObject}
\alias{summary,SaemixObject}
\alias{summary,SaemixObject-method}
\alias{coef,SaemixObject}
\alias{coef,SaemixObject-method}

\alias{SaemixRes-class}
\alias{SaemixRes}
\alias{[<-,SaemixRes-method}
\alias{[,SaemixRes-method}
\alias{initialize,SaemixRes-method}
\alias{print,SaemixRes-method}
\alias{showall,SaemixRes-method}
\alias{show,SaemixRes-method}
\alias{print,SaemixRes}
\alias{showall,SaemixRes}
\alias{show,SaemixRes}

\title{Class "SaemixObject"}
\description{ An object of the SaemixObject class, storing the results obtained by a call to the SAEM algorithm
}
\section{Objects from the Class}{
Objects are created by a call \code{saemix()}.
}
\section{Slots}{
  \describe{
    \item{\code{data}:}{Object of class \code{"SaemixData"} an object of the SaemixData class, representing a longitudinal data structure, used by the SAEM algorithm. See \code{\link{SaemixData}}}
    \item{\code{model}:}{Object of class \code{"SaemixModel"} an object of the SaemixModel class, representing the structure a non-linear mixed effect model, used by the SAEM algorithm. See \code{\link{SaemixModel}} }
    \item{\code{results}:}{Object of class \code{"SaemixRes"} ~~ }
    \item{\code{rep.data}:}{Object of class \code{"SaemixRepData"} (internal use only) the data part, replicated a number of times equal to the number of chains used in the SAEM algorithm (see documentation for details). Not intended to be accessed to directly by the user.}
    \item{\code{sim.data}:}{Object of class \code{"SaemixSimData"} (internal use only) data simulated according to the design in data, with the model in model and the parameters estimated by the SAEM algorithm, after a call to simul.saemix (see documentation for details). Not intended to be accessed to directly by the user. }
    \item{\code{options}:}{Object of class \code{"list"} a list of options containing variables controlling the algorithm }
    \item{\code{prefs}:}{Object of class \code{"list"} a list of graphical preferences applied to plots }
  }
}
\section{Methods}{
  \describe{
    \item{[<-}{\code{signature(x = "SaemixObject")}: replace elements of object }
    \item{[}{\code{signature(x = "SaemixObject")}: access elements of object }
    \item{initialize}{\code{signature(.Object = "SaemixObject")}: internal function to initialise object, not to be used }
    \item{plot}{\code{signature(x = "SaemixObject")}: plot results (see \code{\link{saemix.plot.data}} }
    \item{predict}{\code{signature(object = "SaemixObject")}: compute model predictions }
    \item{print}{\code{signature(x = "SaemixObject")}: prints details about the object }
    \item{showall}{\code{signature(object = "SaemixObject")}: prints an extensive summary of the object }
    \item{show}{\code{signature(object = "SaemixObject")}: prints a short summary of the object }
    \item{summary}{\code{signature(object = "SaemixObject")}: summary of the results}
    \item{coef}{\code{signature(object = "SaemixObject")}: extracts coefficients. Returns a list with components fixed (estimated fixed effects), population (population parameter estimates, including covariate effects: a list with two components map and cond), individual (individual parameter estimates: a list with two components map and cond). For population and individual, the map component of the list gives the MAP estimates (the mode of the distribution) while the cond component gives the conditional mean estimates. Some components may be missing (eg MAP estimates) if they have not been computed during or after the fit.}
	 }
}
\details{ Details of the algorithm can be found in the pdf file accompanying the package.}
\references{ Kuhn E, Lavielle M. Maximum likelihood estimation in nonlinear mixed effects models. Computational Statistics and Data Analysis 49, 4 (2005), 1020-1038.

Comets E, Lavenu A, Lavielle M. SAEMIX, an R version of the SAEM algorithm. 20th meeting of the Population Approach Group in Europe, Athens, Greece (2011), Abstr 2173.
}
\author{ Emmanuelle Comets <emmanuelle.comets@inserm.fr>, Audrey Lavenu, Marc Lavielle.
}
\seealso{\code{\link{SaemixData}},\code{\link{SaemixModel}}, \code{\link{saemixControl}},\code{\link{saemix}},\code{\link{plot.saemix}}, \code{\link{saemix.plot.data}}}
\examples{
showClass("SaemixObject")
}
\keyword{classes}
