% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{transformDF}
\alias{transformDF}
\title{Transform Data Frame}
\usage{
transformDF(df)
}
\arguments{
\item{df}{The input data frame in the usual PanelPRO format.}
}
\value{
A data frame in the format required for clipp with the following columns:
\item{individual}{ID of the individual}
\item{isProband}{Indicator if the individual is a proband}
\item{family}{Family ID}
\item{mother}{Mother's ID}
\item{father}{Father's ID}
\item{aff}{Affection status}
\item{sex}{Sex (2 for female, 1 for male)}
\item{age}{Age at diagnosis or current age}
\item{geno}{Genotype information}
}
\description{
This function transforms a data frame from the standard format used in PanelPRO
into the required format which conforms to the requirements of penetrance (and clipp).
}
\examples{
# Create example data frame
df <- data.frame(
  ID = 1:2,
  PedigreeID = c(1,1),
  Sex = c(0,1),
  MotherID = c(NA,1),
  FatherID = c(NA,NA),
  isProband = c(1,0),
  CurAge = c(45,20),
  isAff = c(1,0),
  Age = c(40,NA),
  Geno = c(1,0)
)

# Transform the data frame
transformed_df <- transformDF(df)
}
