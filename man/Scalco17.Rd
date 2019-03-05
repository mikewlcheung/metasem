\name{Scalco17}
\alias{Scalco17}
\docType{data}
\title{Correlation Matrices from Scalco et al. (2017)}
\description{
  The data set includes correlation matrices using the theory of planned
  behavior to predict organic food consumption reported by Scalco17 et al. (2017). 
}
\usage{data(Scalco17)}

\details{
A list of data with the following structure:
\describe{
  \item{data}{A list of correlation matrices. The variables are
	\emph{ATT} (attitude), \emph{SN} (subjective norm), \emph{PBC}
	(perceived behavior control), \emph{BI} (behavioral intention), and \emph{BEH} (behavior)}
  \item{n}{A vector of sample sizes}
  \item{Age}{A vector of the mean age of the samples}
  \item{Female}{A vector of the percentage of the female samples}
}
}
\source{
  Scalco, A., Noventa, S., Sartori, R., & Ceschi, A. (2017). Predicting organic food consumption: A meta-analytic structural equation model based on the theory of planned behavior. \emph{Appetite}, \bold{112}, 235-248.
}

\examples{
data(Scalco17)
}
\keyword{datasets}
