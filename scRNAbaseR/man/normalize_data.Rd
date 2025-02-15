% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scripts.R
\name{normalize_data}
\alias{normalize_data}
\title{Normalize gene expression data}
\usage{
normalize_data(data, type = "counts", factor = 1e+06)
}
\arguments{
\item{data}{A matrix or data frame of raw gene expression data.}

\item{type}{A type of normalization based on number of positive expressed genes or total counts per cell / set. Default is 'counts'}

\item{factor}{A numeric scaling factor for normalization. Default is 1,000,000.}
}
\value{
A normalized dataset with log2-transformed values.
}
\description{
Applies log2 normalization to gene expression data based on gene counts and a scaling factor.
}
\examples{
normalized_data <- normalize_data(data)
}
