\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot,ANY-method}
\title{ Methods for Function plot }
\description{ Methods for function \code{plot}}
\section{Methods}{
\describe{

\item{\code{signature(x = "ANY")}}{ default plot function ?}

\item{\code{signature(x = "SaemixData")}}{ Plots the data. Defaults to a spaghetti plot of response versus predictor, with lines joining the data for one individual.}

\item{\code{signature(x = "SaemixModel")}}{ Plots prediction of the model }

\item{\code{signature(x = "SaemixObject")}}{ This method gives access to a number of plots that can be performed on a SaemixObject}

\item{\code{signature(x = "SaemixSimData")}}{ Plots simulated datasets}
}}
\keyword{methods}

