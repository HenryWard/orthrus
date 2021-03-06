% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load.R
\name{split_guides_by_type}
\alias{split_guides_by_type}
\title{Splits guides by type.}
\usage{
split_guides_by_type(guides)
}
\arguments{
\item{guides}{Dataframe returned from \code{retrieve_guides_by_gene}.}
}
\value{
A list of three separate guide lists. Combinatorial-targeting guides are
  stored in the key "combn", single-targeting guides are stored in the 
  key "single" and dual-targeting guides are stored in the key "dual".
}
\description{
Splits guide returned from \code{retrieve_guides_by_gene} by targeting type, which is 
one of dual-gene exonic-exonic, single-gene exonic-exonic or exonic-intergenic.
}
