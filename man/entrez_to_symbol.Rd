% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/entrez_to_symbol.R
\name{entrez_to_symbol}
\alias{entrez_to_symbol}
\title{Obtain gene symbols for entrezgene IDs}
\usage{
entrez_to_symbol(df, species, dir_annotation = getwd(), verbose = TRUE)
}
\arguments{
\item{df}{A data.frame containing a column named "entrezgene_id".}

\item{species}{One of "mouse" or "human".}

\item{dir_annotation}{The directory to store annotation reference.}

\item{verbose}{Set to FALSE to avoid messages.}
}
\value{
The data.frame df with an added column containing gene symbols.
}
\description{
Obtain gene symbols for entrezgene IDs
}
