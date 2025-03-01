% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/F_transformation.R
\name{F_transformation}
\alias{F_transformation}
\title{Transformation Function for Mixed Graphical Models}
\usage{
F_transformation(X, type, parameter, cont = FALSE)
}
\arguments{
\item{X}{A numeric matrix or vector.}

\item{type}{A character vector specifying the type of each variable.
Options: "c" (continuous), "d" (discrete), "z" (zero-inflated), "m" (categorical).}

\item{parameter}{Numeric parameter used in the transformation.}

\item{cont}{Logical. If TRUE, continuous variables are transformed. Default is FALSE.}
}
\value{
A matrix of transformed data.
}
\description{
Performs transformation of input data based on its type.
Used for continuous, discrete, and zero-inflated data in Mixed Graphical Models.
}
