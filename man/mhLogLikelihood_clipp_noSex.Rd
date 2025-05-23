% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mhLoglikehood.r
\name{mhLogLikelihood_clipp_noSex}
\alias{mhLogLikelihood_clipp_noSex}
\title{Calculate Log Likelihood without Sex Differentiation}
\usage{
mhLogLikelihood_clipp_noSex(
  paras,
  families,
  twins,
  max_age,
  baseline_data,
  prev,
  geno_freq,
  trans,
  BaselineNC,
  ncores
)
}
\arguments{
\item{paras}{Numeric vector, the parameters for the Weibull distribution and scaling factors.
Should contain in order: gamma, delta, given_median, given_first_quartile.}

\item{families}{Data frame, containing pedigree information with columns for 'age', 'aff' (affection status), and 'geno' (genotype).}

\item{twins}{Information on monozygous twins or triplets in the pedigrees.}

\item{max_age}{Integer, maximum age considered in the analysis.}

\item{baseline_data}{Numeric vector, baseline risk data for each age.}

\item{prev}{Numeric, prevalence of the risk allele in the population.}

\item{geno_freq}{Numeric vector, represents the frequency of the risk type and its complement in the population.}

\item{trans}{Numeric matrix, transition matrix that defines the probabilities of allele transmission from parents to offspring.}

\item{BaselineNC}{Logical, indicates if non-carrier penetrance should be based on the baseline data or the calculated non-carrier penetrance.}

\item{ncores}{Integer, number of cores to use for parallel computation.}
}
\value{
Numeric, the calculated log likelihood.
}
\description{
This function calculates the log likelihood for a set of parameters and data without considering sex differentiation using the clipp package.
}
\references{
Details about the clipp package and methods can be found in the package documentation.
}
