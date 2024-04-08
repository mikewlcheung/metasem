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
              cor.analysis = TRUE, acov="weighted", append.vars=TRUE,
              asyCovOld=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A list of data with correlation/covariance matrix in \code{x$data} and
   sample sizes \code{x$n}. Additional variables in \code{x} can be attached.
}
\item{n}{If \code{x} is a list of correlation matrices without
  \code{x$data} and \code{x$n}, a vector of sample sizes \code{n} must
  be provided.}
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
\item{append.vars}{Whether to append the additional variables to
  the output dataframe.}
\item{asyCovOld}{Whether to use the old version of \code{asyCov}. See \code{\link[metaSEM]{asyCov}}.}
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
\donttest{
## Provide a list of correlation matrices and a vector of sample sizes as the inputs
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
my.df2 <- Cor2DataFrame(Nohe15A1)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ osmasem }
