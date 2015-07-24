\name{SaemixModel-class}
\Rdversion{1.1}
\docType{class}
\alias{SaemixModel-class}
\alias{SaemixModel}
\alias{[<-,SaemixModel-method}
\alias{[,SaemixModel-method}
\alias{initialize,SaemixModel-method}
\alias{plot,SaemixModel-method}
\alias{print,SaemixModel-method}
\alias{showall,SaemixModel-method}
\alias{show,SaemixModel-method}
\alias{plot,SaemixModel}
\alias{print,SaemixModel}
\alias{showall,SaemixModel}
\alias{show,SaemixModel}
\alias{summary,SaemixModel}
\alias{summary,SaemixModel-method}

\title{Class "SaemixModel"}
\description{ An object of the SaemixModel class, representing a non-linear mixed-effect model structure, used by the SAEM algorithm.
}
\section{Objects from the Class}{An object of the SaemixModel class can be created by using the function \code{\link{saemixModel}}
}
\section{Slots}{
  \describe{
    \item{\code{model}:}{ function used to compute the structural model. The function should return a vector of predicted values given a matrix of individual parameters, a vector of indices specifying which records belong to a given individual, and a matrix of dependent variables (see examples).}
    \item{\code{description}:}{ model description (optional) as a character string}
    \item{\code{psi0}:}{ a matrix with a number of columns equal to the number of parameters in the model, and one (when no covariates are available) or two (when covariates enter the model) giving the initial estimates for the fixed effects. The column names of the matrix should be the names of the parameters in the model, and will be used in the plots and the summaries}
    \item{\code{transform.par}:}{ the distribution for each parameter (0=normal, 1=log-normal, 2=probit, 3=logit). Defaults to a vector of 1s (all parameters have a log-normal distribution)}
    \item{\code{fixed.estim}:}{ whether parameters should be estimated (1) or fixed to their initial estimate (0). Defaults to a vector of 1s}
    \item{\code{error.model}:}{ name of the residual error model}
    \item{\code{covariate.model}:}{ a matrix giving the covariate model. Defaults to no covariate in the model (empty matrix)}
    \item{\code{betaest.model}:}{ a matrix giving the effects model (internal)}
    \item{\code{covariance.model}:}{ a square matrix of size equal to the number of parameters in the model, giving the variance-covariance matrix of the model: 1s correspond to estimated variances (in the diagonal) or covariances (off-diagonal elements). Defaults to the identity matrix}
    \item{\code{omega.init}:}{ a square matrix of size equal to the number of parameters in the model, giving the initial estimate for the variance-covariance matrix of the model. Defaults to the identity matrix}
    \item{\code{error.init}:}{ a vector of size 2 giving the initial value of a and b in the error model. Defaults to 1 for each estimated parameter in the error model}
    \item{\code{nb.parameters}:}{ number of parameters}
    \item{\code{name.modpar}:}{ names of the model parameters}
    \item{\code{name.fixed}:}{ names of the fixed effects estimated in the model}
    \item{\code{name.random}:}{ names of the random effects estimated in the model}
    \item{\code{name.res}:}{ names of the parameters of the residual error model}
    \item{\code{name.predictors}:}{ names of the predictors (X)}
    \item{\code{name.X}:}{ name of the predictor used in graphs}
    \item{\code{name.response}:}{ name of the response (Y)}
    \item{\code{name.cov}:}{ name of the covariates}
    \item{\code{indx.fix}:}{ index of estimated fixed effects (internal)}
    \item{\code{indx.cov}:}{ index of estimated fixed effects associated with covariate effects (internal)}
    \item{\code{indx.omega}:}{ index of estimated random effects (internal)}
    \item{\code{indx.res}:}{ index of parameters of the residual error model (internal)}
    \item{\code{Mcovariates}:}{ matrix of the covariates (internal)}
  }
}
\section{Methods}{
  \describe{
    \item{[<-}{\code{signature(x = "SaemixModel")}: replace elements of object }
    \item{[}{\code{signature(x = "SaemixModel")}: access elements of object }
    \item{initialize}{\code{signature(.Object = "SaemixModel")}: internal function to initialise object, not to be used }
    \item{plot}{\code{signature(x = "SaemixModel")}: plot results (see \code{\link{saemix.plot.data}} }
    \item{print}{\code{signature(x = "SaemixModel")}: prints details about the object }
    \item{showall}{\code{signature(object = "SaemixModel")}: prints an extensive summary of the object }
    \item{show}{\code{signature(object = "SaemixModel")}: prints a short summary of the object }
    \item{summary}{\code{signature(object = "SaemixModel")}: summary of the model}
	 }
}
\references{ Kuhn E, Lavielle M. Maximum likelihood estimation in nonlinear mixed effects models. Computational Statistics and Data Analysis 49, 4 (2005), 1020-1038.

Comets E, Lavenu A, Lavielle M. SAEMIX, an R version of the SAEM algorithm. 20th meeting of the Population Approach Group in Europe, Athens, Greece (2011), Abstr 2173.
}
\author{ Emmanuelle Comets <emmanuelle.comets@inserm.fr>, Audrey Lavenu, Marc Lavielle.
}
\seealso{\code{\link{SaemixData}},\code{\link{saemixModel}}, \code{\link{saemixControl}},\code{\link{saemix}}}
\examples{
showClass("SaemixModel")
}
\keyword{classes}
