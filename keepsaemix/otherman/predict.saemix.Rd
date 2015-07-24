\name{predict.saemix}
\alias{predict.saemix}
\title{ SAEM model predictions
}
\description{
Obtain model predictions using the results of an SAEM fit
}
\usage{
predict.saemix(saemix.res,saemix.data,saemix.model,saemix.options)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{saemix.res}{ saemix.res: corresponding element of the list returned by fit.saemix, containing the model results}
  \item{saemix.data}{ saemix.data: corresponding element of the list returned by fit.saemix, containing the data file and identifiers}
  \item{saemix.model}{ saemix.model: corresponding element of the list returned by fit.saemix, containing the characteristics of the estimated model and parameters}
  \item{saemix.options}{ saemix.options: corresponding element of the list returned by fit.saemix, containing the fit options}
}
\details{
A call to this function is required before requesting results to be plotted. If \code{\link{map.saemix}} has not yet been run, it will be called to compute the individual parameter estimates necessary to obtain model predictions. Predictions are obtained for both MAP and conditional estimates.
}
\value{ Returns an updated list saemix.res, including components
\describe{
  \item{ipred}{model predictions obtained using the individual MAP parameters (map)}
  \item{ipred.cond}{model predictions obtained using the conditional mean parameters (cond.mean.phi)}
  \item{ppred}{The model predictions obtained using the population parameters ( f(x_i,h(C_i,hat{mu})) )}
  \item{iwres}{individual standardised residuals, using the model predictions obtained with the MAP individual estimates of the parameters}
  \item{icwres}{individual standardised residuals, using the model predictions obtained with the conditional individual estimates of the parameters}
}
}
\references{ Kuhn, E., and Lavielle, M. Maximum likelihood estimation in nonlinear mixed effects models. Computational Statistics and Data Analysis 49, 4 (2005), 1020-1038.

Monolix32_UsersGuide.pdf (\url{http://software.monolix.org/sdoms/software/})
}
\author{ Emmanuelle Comets <emmanuelle.comets@inserm.fr>, Audrey Lavenu, Marc Lavielle.
}
\seealso{\code{\link{SaemixObject}},\code{\link{saemix}}}, \code{\link{saemix.plot.data}} }
\examples{

# Not run: update saemix.res after an SAEM fit (see eg saemix.plot.data)
# saemix.res<-predict.saemix(saemix.res,saemix.data,saemix.model,saemix.options)

}
\keyword{model}
