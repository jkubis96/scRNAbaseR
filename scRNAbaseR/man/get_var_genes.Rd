% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scripts.R
\name{get_var_genes}
\alias{get_var_genes}
\title{Extract variable or equally expressed genes}
\usage{
get_var_genes(var_data, side = "equal")
}
\arguments{
\item{var_data}{A data frame of gene variance and mean values.}

\item{side}{A string specifying which genes to extract ('equal' or 'variable'). Default is 'equal'.}
}
\value{
A data frame of selected genes based on the specified criterion.
}
\description{
Selects genes based on their variance and mean expression values.
}
\examples{
variable_genes <- get_var_genes(var_data, side = 'variable')
}
