\name{Cor2DataFrame}
\alias{Cor2DataFrame}
\alias{tssem2DataFrame}
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
tssem2DataFrame(X, v.na.replace = TRUE, row.names.unique = FALSE,
                cor.analysis = TRUE, acov="weighted", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A correlation/covariance matrix or a list of
	correlation/covariance matrices.}
   \item{X}{A list of correlation/covariance matrix in \code{X$data} and
   a list of sample sizes \code{X$n}. More variables can be included in
   \code{X}. These variables will be added in the output's \code{data}.} 
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

\value{A list of components: (1) a data frame of correlations or covariances with their
  sampling covariance matrices; (2) a vector of sample
  sizes; (3) labels of the correlations; and (3) labels of their sampling covariance matrices. }

\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\seealso{ \code{\link[metaSEM]{asyCov}}, \code{\link[metaSEM]{osmasem}}, \code{\link[metaSEM]{create.vechsR}},
  \code{\link[metaSEM]{create.Tau2}}, \code{\link[metaSEM]{create.V}}
}


\examples{
\dontrun{   
my.df1 <- Cor2DataFrame(Nohe15A1$data, Nohe15A1$n)

## Add Lag time as a variable
my.df1$data <- data.frame(my.df1$data, Lag=Nohe15A1$Lag, check.names=FALSE)

## Data
my.df1$data

## Sample sizes
my.df1$n

## ylabels
my.df1$ylabels

## vlabels
my.df1$vlabels

#### Simplified version to do it
my.df2 <- tssem2DataFrame(Nohe15A1)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ osmasem }
