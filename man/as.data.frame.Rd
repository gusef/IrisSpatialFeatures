% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transform.functions.R
\docType{methods}
\name{as.data.frame,ImageSet-method}
\alias{as.data.frame,ImageSet-method}
\title{Read inForm output and store it in an IrisSpatialFeatures ImageSet object.}
\usage{
\S4method{as.data.frame}{ImageSet}(x)
}
\arguments{
\item{x}{Iris ImageSet boject}

\item{...}{Additional arguments}
}
\value{
A dataframe
}
\description{
Read inForm output and store it in an IrisSpatialFeatures ImageSet object.
}
\examples{
dataset <- IrisSpatialFeatures_data
dataframe <- as.data.frame(dataset)
}