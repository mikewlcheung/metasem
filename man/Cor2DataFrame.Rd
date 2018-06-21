\name{Cor2DataFrame}
\alias{Cor2DataFrame}
\title{Convert correlation or covariance matrices into a dataframe of correlations or
  covariances with their sampling covariance matrices
}
\description{It converts the correlation or covariance matrices into a
  dataframe of correlations or covariances with their asymptotic
  sampling covariance matrices. It uses the \code{asyCov} at the backend. 
}
\usage{
Cor2DataFrame(x, n, v.na.replace = TRUE, row.names.unique = FALSE,
              cor.analysis = TRUE, acov="weighted", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A correlation/covariance matrix or a list of
	correlation/covariance matrices.}
\item{n}{Sample size or a vector of sample sizes}
\item{v.na.replace}{Logical. Missing value is not allowed in definition
  variables. If it is \code{TRUE} (the default), missing value is
  replaced by a large value (1e10). These values are not used in the analysis.} 
\item{row.names.unique}{Logical, If it is \code{FALSE} (the default), unique
  row names are not created.}
\item{cor.analysis}{Logical. The output is either a correlation or
  covariance matrix.}
\item{acov}{If it is \code{weighted}, the average correlation/covariance
  matrix is calculated based on the weighted mean with the sample
  sizes. The average correlation/covariance matrix is used to calculate the sampling
	variance-covariance matrices.}
\item{\dots}{Further arguments to be passed to \code{\link[metaSEM]{asyCov}}.}
}

\value{A list of components: (1) a dataframe of correlations or covariances with their
  sampling covariance matrices; (2) a vector of sample
  sizes; (3) labels of the correlations; and (3) lables of their sampling covariance matrices. }

\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\seealso{ \code{\link[metaSEM]{asyCov}}, \code{\link[metaSEM]{osmasem}}, \code{\link[metaSEM]{create.vechsR}},
  \code{\link[metaSEM]{create.Tau2}}, \code{\link[metaSEM]{create.V}}
}


\examples{
\dontrun{
my.df <- Cor2DataFrame(Nohe15A1$data, Nohe15A1$n)

## Data
my.df$data

## Sample sizes
my.df$n

## ylabels
my.df$ylabels

## vlabels
my.df$vlabels
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ osmasem }
