% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mhLoglikehood.r
\name{lik_noSex}
\alias{lik_noSex}
\title{Likelihood Calculation without Sex Differentiation}
\usage{
lik_noSex(
  i,
  data,
  alpha,
  beta,
  delta,
  gamma,
  max_age,
  baselineRisk,
  BaselineNC,
  prev
)
}
\arguments{
\item{i}{Integer, index of the individual in the data set.}

\item{data}{Data frame, containing individual demographic and genetic information. Must include columns for 'age', 'aff' (affection status), and 'geno' (genotype).}

\item{alpha}{Numeric, Weibull distribution shape parameter.}

\item{beta}{Numeric, Weibull distribution scale parameter.}

\item{delta}{Numeric, shift parameter for the Weibull function.}

\item{gamma}{Numeric, asymptote parameter (only scales the entire distribution).}

\item{max_age}{Integer, maximum age considered in the analysis.}

\item{baselineRisk}{Numeric vector, baseline risk for each age.}

\item{BaselineNC}{Logical, indicates if non-carrier penetrance should be based on SEER data or the calculated non-carrier penetrance.}

\item{prev}{Numeric, prevalence of the risk allele in the population.}
}
\value{
Numeric vector, containing likelihood values for unaffected and affected individuals.
}
\description{
This function calculates the likelihood for an individual based on Weibull distribution parameters without considering sex differentiation.
}
