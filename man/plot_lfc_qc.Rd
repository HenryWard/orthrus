% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_lfc_qc}
\alias{plot_lfc_qc}
\title{Plot replicate comparisons.}
\usage{
plot_lfc_qc(
  df,
  screens,
  output_folder,
  negative_controls = NULL,
  control_label = "Negative controls",
  plot_type = "png",
  display_numbers = TRUE
)
}
\arguments{
\item{df}{Reads or lfc dataframe.}

\item{screens}{List of screens created with \code{add_screens}.}

\item{output_folder}{Folder to output plots to.}

\item{negative_controls}{List of negative control genes to append to default list of
non-essential genes (default NULL).}

\item{control_label}{Label for negative control genes in plot (default "Negative controls").}

\item{plot_type}{Type of plot to output, one of "png" or "pdf" (default "png").}

\item{display_numbers}{Whether or not to include PCC values in heatmap (default TRUE).}
}
\description{
Plots replicate comparisons for all replicates in a list of screens and outputs
plots to a given folder.
}
