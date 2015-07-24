\name{print-methods}
\docType{methods}
\alias{print.saemix}
\alias{print-methods}
\alias{print,ANY-method}
\title{ Methods for Function print }
\description{ Prints a summary of an object}
\section{Methods}{
\describe{

\item{\code{signature(x = "ANY")}}{
Default print function
}

\item{\code{signature(x = "SaemixData")}}{
Prints a summary of a SaemixData object
}

\item{\code{signature(x = "SaemixModel")}}{
Prints a summary of a SaemixModel object
}

\item{\code{signature(x = "SaemixObject")}}{
Prints a summary of the results from a SAEMIX fit
}

\item{\code{signature(x = "SaemixRes")}}{
Not user-level
}
}}
\keyword{methods}
