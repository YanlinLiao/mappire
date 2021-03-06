% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aiding_functions.R
\name{found_probability}
\alias{found_probability}
\title{found_probability
Find the possible gamete genotypes which contribute by each parents and their corresponding
probabilities.}
\usage{
found_probability(
  P1,
  P2,
  inheritance1,
  inheritance2,
  row = F1[i, ],
  ind_nme = F1_ind[i],
  ploidy
)
}
\arguments{
\item{P1}{parent1's genotype}

\item{P2}{parent2's genotype}

\item{inheritance1}{inheritance of P1, can be either 'Polysomic' or 'Disomic'}

\item{inheritance2}{inheritance of P2, can be either 'Polysomic' or 'Disomic'}

\item{row}{genotype from one offspring, eg. "A" "A" "E" "F"}

\item{ind_nme}{the name of that offspring, eg. "F1_4"}
}
\value{
a dataframe contains the possible offspring genotype and their corresponding
probability, the colnames including:
'individual P1_1 P1_2 P2_1 P2_2 probability   pairing pair_freq P1_gamete_freq P2_gamete_freq '
    ' F1_4    A    A    E    F           1 AAAB_EEEF         1            0.5            0.5 '
}
\description{
found_probability
Find the possible gamete genotypes which contribute by each parents and their corresponding
probabilities.
}
