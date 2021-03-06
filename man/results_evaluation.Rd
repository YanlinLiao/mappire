% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exported_functions.R
\name{results_evaluation}
\alias{results_evaluation}
\title{evaluate estimated linkage result}
\usage{
results_evaluation(
  result = Linkage_res$output,
  genfile = genfile,
  mapfile = mapfile,
  plot = TRUE,
  filename
)
}
\arguments{
\item{result}{linkage estimation results. A dataframe include:
marker_P, marker_Q, r, LOD, phasing_P, phasing_Q}

\item{genfile}{simulated .gen file used for check}

\item{mapfile}{simulated .map file used for check}

\item{plot}{TRUE/FALSE. If True, a plot will be saved in your folder 'Plot/'}

\item{filename}{filename that the figure will be saved}
}
\value{
a dataframe include the linkage results and simulated r and phasing
}
\description{
evaluate estimated linkage result
}
