% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/counting.functions.R
\docType{methods}
\name{counts_data_frame}
\alias{counts_data_frame}
\alias{counts_data_frame,ImageSet-method}
\alias{counts_data_frame,ANY,ANY-method}
\title{DataFrame from all the counts per frame}
\usage{
counts_data_frame(x, ...)

\S4method{counts_data_frame}{ImageSet}(x)
}
\arguments{
\item{x}{IrisSpatialFeatures ImageSet object.}

\item{...}{Additional arguments}
}
\value{
data frame
}
\description{
DataFrame from all the counts per frame
}
\examples{

#loading pre-read dataset
dataset <- IrisSpatialFeatures_data
counts_data_frame(dataset)

}
