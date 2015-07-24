\name{saemixData}
\alias{saemixData}
\title{ Function to create a SaemixData object
}
\description{ This function creates a SaemixData object. The only mandatory argument is the name of the dataset. If the dataset has a header (or named columns), the program will attempt to detect which column correspond to ID, predictor(s) and response. Warning messages will be printed during the object creation and should be read for details.
}
\usage{
saemixData(name.data, header, sep, na, name.group, name.predictors,
  name.response, name.X, name.covariates = c(),
  units = list(x = "", y = "", covariates = c()))
}
\arguments{
    \item{name.data}{ name of the dataset (can be a character string giving the name of a file on disk or of a dataset in the R session, or the name of a dataset}
    \item{header}{ whether the dataset/file contains a header. Defaults to TRUE }
    \item{sep}{ the field separator character. Defaults to any number of blank spaces ("") }
    \item{na}{ a character vector of the strings which are to be interpreted as NA values. Defaults to c(NA) }
    \item{name.group}{ name (or number) of the column containing the subject id}
    \item{name.predictors}{ name (or number) of the column(s) containing the predictors (the algorithm requires at least one predictor x)}
    \item{name.response}{ name (or number) of the column containing the response variable y modelled by predictor(s) x}
    \item{name.covariates}{ name (or number) of the column(s) containing the covariates, if present (otherwise missing)}
    \item{name.X}{ name of the column containing the regression variable to be used on the X axis in the plots (defaults to the first predictor)}
    \item{units}{ list with up to three elements, x, y and optionally covariates, containing the units for the X and Y variables respectively, as well as the units for the different covariates (defaults to empty)}
}
\details{ This function is the user-friendly constructor for the SaemixData object class.
}
\value{ A SaemixData object (see \code{\link{saemixData}}).
}
\references{ Kuhn, E., and Lavielle, M. Maximum likelihood estimation in nonlinear mixed effects models. Computational Statistics and Data Analysis 49, 4 (2005), 1020-1038.

Monolix32_UsersGuide.pdf (\url{http://software.monolix.org/sdoms/software/})
}
\author{ Emmanuelle Comets <emmanuelle.comets@inserm.fr>, Audrey Lavenu, Marc Lavielle.
}
\seealso{\code{\link{SaemixData}},\code{\link{SaemixModel}}, \code{\link{saemixControl}},\code{\link{saemix}}}
\examples{

data(theo.saemix)

saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, 
  name.group=c("Id"),name.predictors=c("Dose","Time"),
  name.response=c("Concentration"),name.covariates=c("Weight","Sex"),
  units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")

print(saemix.data)

plot(saemix.data)
}
\keyword{models}
