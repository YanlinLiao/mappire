% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exported_functions.R
\name{all_equations_generator}
\alias{all_equations_generator}
\title{Generate all possible equations that will be used in your datasets
Here it made use of the generalized equation and your own dataset to create all possible equations that probably used later.}
\usage{
all_equations_generator(
  P1_inheritance = "Polysomic",
  P2_inheritance = "Disomic",
  MT = MT,
  ploidy
)
}
\arguments{
\item{P1_inheritance}{'polysomic' or 'disomic'. It will depend on this to generate different equations}

\item{P2_inheritance}{'polysomic' or 'disomic'. It will depend on this to generate different equations}

\item{MT}{a dataframe include markername, parent genotype, P1 and P2 marker type, the number of
shared haplotypes, the number of unique haplotypes, combined marker type}

\item{ploidy}{crop ploidy level}
}
\value{
a list include:
 - equations: equations for P1 and P2 separately
 - input: P1 and P2's inheritance, ploidy level, and marker type (MT)
}
\description{
Generate all possible equations that will be used in your datasets
Here it made use of the generalized equation and your own dataset to create all possible equations that probably used later.
}
