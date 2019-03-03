\name{create.mxModel}
\alias{create.mxModel}
\title{Creating mxModels
}
\description{It creates an mxModel from an RAM object.  
}
\usage{
create.mxModel(model.name="mxModel", RAM=NULL, Amatrix=NULL,
               Smatrix=NULL, Fmatrix=NULL, Mmatrix=NULL,
               data, intervals.type = c("z", "LB"),
               mx.algebras=NULL, mxModel.Args=NULL,
               mxRun.Args=NULL, suppressWarnings=TRUE,
               silent=TRUE, run=TRUE, ...)
}
\arguments{
  \item{model.name}{A string for the model name in \code{\link[OpenMx]{mxModel}}. 
}
\item{RAM}{An RAM object including a list of matrices of the model
  returned from \code{\link[metaSEM]{lavaan2RAM}}.}
\item{Amatrix}{If \code{RAM} is not specified, an \code{Amatrix} is
  required. Please note that the dimension names are required in the
  \code{Amatrix$labels} when \code{RAM} is not included.}
\item{Smatrix}{If \code{RAM} is not specified, an \code{Smatrix} is
  required.}
\item{Fmatrix}{If \code{RAM} is not specified, an \code{Fmatrix} is
  required.}
\item{Mmatrix}{If \code{RAM} is not specified, an \code{Mmatrix} is
  required.}
\item{data}{A data frame for the data.}
\item{intervals.type}{Either \code{z} (default if missing) or
  \code{LB}. If it is \code{z}, it calculates the 95\% confidence
  intervals (CIs) based on the estimated standard error. If it
  is \code{LB}, it calculates the 95\% likelihood-based CIs on the parameter estimates.}
\item{mxModel.Args}{A list of arguments passed to \code{\link[OpenMx]{mxModel}}.}
\item{mxRun.Args}{A list of arguments passed to \code{\link[OpenMx]{mxRun}}.}
\item{suppressWarnings}{Logical. If it is \code{TRUE}, warnings are
	suppressed. This argument is passed to \code{\link[OpenMx]{mxRun}}.}
\item{silent}{Logical. The argument is passed to \code{\link[OpenMx]{mxRun}}}
\item{run}{Logical. If \code{FALSE}, only return the mx model without running the analysis.}
\item{\dots}{Not used yet.}  
}

\value{An object of class \code{mxModel}}

\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}

\examples{
\dontrun{
## Generate data
set.seed(100)
n <- 100
x <- rnorm(n)
y <- 0.5*x + rnorm(n, mean=0, sd=sqrt(1-0.5^2))
my.df <- data.frame(y=y, x=x)

## Regression model
model <- "y ~ x   # Regress y on x
          y ~ 1   # Intercept of y
          x ~ 1   # Mean of x"
plot(model)

RAM <- lavaan2RAM(model, obs.variables=c("y", "x"))

my.fit <- create.mxModel(RAM=RAM, data=my.df)
summary(my.fit)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
