% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/score.R
\name{score_combn_vs_single}
\alias{score_combn_vs_single}
\title{Scores conditions against a single control.}
\usage{
score_combn_vs_single(
  combn_guides,
  single_guides,
  screens,
  screen_names,
  min_guides = 3,
  test = "moderated-t",
  loess = TRUE,
  fdr_method = "BY",
  filter_genes = NULL,
  ignore_orientation = FALSE,
  collapse_single_targeting = TRUE,
  return_residuals = TRUE,
  verbose = FALSE
)
}
\arguments{
\item{combn_guides}{A list of exonic-exonic guides returned from \code{split_guides}.}

\item{single_guides}{A list of exonic-intergenic guides returned from \code{split_guides}.}

\item{screens}{List of screens generated with \code{add_screens}.}

\item{screen_names}{A list of screen names to score against a derived null model from
single-gene effects.}

\item{min_guides}{The minimum number of guides per gene pair required to score data 
(default 3).}

\item{test}{Type of hypothesis testing to run. Must be one of "rank-sum" for Wilcoxon
rank-sum testing or "moderated-t" for moderated t-testing (default "moderated-t").}

\item{loess}{If true, loess-normalizes residuals before running hypothesis testing.
Only works when test = "moderated-t" (default TRUE).}

\item{fdr_method}{Type of FDR to compute. One of "BH", "BY" or "bonferroni" (default "BY").}

\item{filter_genes}{List of genes to filter from scoring (default NULL).}

\item{ignore_orientation}{If TRUE, aggregates guides across both orientations, returning only
one p-value and FDR column with orientation2 p-values set to NA (default FALSE).}

\item{collapse_single_targeting}{If TRUE, takes the mean of single-targeting controls when
there are multiple controls that match a given gene pair (default TRUE).}

\item{return_residuals}{If FALSE, doesn't return residuals dataframe (default TRUE).
This is recommend if scoring large datasets and memory is a limitation.}

\item{verbose}{If true, prints verbose output (default FALSE).}
}
\value{
A list containing two dataframes. The first entry, named "scored_data" in the list,
  contains scored data with separate columns given by the specified control and condition
  names. The second entry, named "residuals" in the list, is a dataframe containing control,
  condition and loess-normalized residuals for all guides.
}
\description{
Scores guides for any number of condition screens against a multiplicative null 
model derived from single-gene effects (e.g. for comparing paralog double-knockouts
to single-knockouts for each paralogous gene). Do NOT use this function if you 
would like to directly compare conditions against a control. In that case, call 
\code{score_conditions_vs_control} instead. After running this function, pass the 
resulting dataframe to \code{call_combn_hits} to call significant effects.
}
