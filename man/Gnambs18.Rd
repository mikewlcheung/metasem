\name{Gnambs18}
\alias{Gnambs18}
\docType{data}
\title{Correlation Matrices from Gnambs, Scharl, and Schroeders (2018)}
\description{
  The data set includes 113 correlation matrices on the Rosenberg
  Self-Esteem Scale reported by Gnambs, Scharl, and Schroeders
  (2018). Thirty-six studies were based on the reported correlation
  matrices (\code{CorMat=1}) whereas the correlation matrices of the
  other 77 studies were calculated from the reported factor loadings.
}
\usage{data(Gnambs18)
}
\details{
A list of data with the following structure:
\describe{
  \item{data}{A list of 113 correlation matrices. The variable names are from \emph{I1} to \emph{I10}.}
  \item{n}{A vector of sample sizes.}
  \item{Year}{The year of publications.}
  \item{Country}{The country of studies conducted.}
  \item{Language}{The language used in the studies.}
  \item{Publication}{Whether the studies were published (1) or
	unpublished (0).}
  \item{MeanAge}{Mean age of the participants.}
  \item{FemaleProp}{Proportion of the female participants.}
  \item{Individualism}{Individualism score of the country.}
  \item{CorMat}{Whether the correlation matrices are obtained from the original
	studies (1) or reproduced from the factor loadings (0).}
}
}
\source{
  Gnambs, T., Scharl, A., & Schroeders, U. (2018). The structure of the Rosenberg Self-Esteem Scale. \emph{Zeitschrift Fur Psychologie}, \bold{226}(1), 14-29. https://doi.org/10.1027/2151-2604/a000317
}
\keyword{datasets}
