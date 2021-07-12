\name{create.mxModel}
\alias{create.mxModel}
\title{Create an mxModel
}
\description{It creates a mxModel from a RAM object.  
}
\usage{
create.mxModel(model.name="mxModel", RAM=NULL, data,
                           intervals.type = c("z", "LB"),
                           var.names=NULL, mxModel.Args=NULL,
                           run=TRUE, mxTryHard=FALSE,
                           silent=TRUE, ...)
}
\arguments{
  \item{model.name}{A string for the model name in \code{\link[OpenMx]{mxModel}}. 
}
\item{RAM}{A RAM object including a list of matrices of the model
  returned from \code{\link[metaSEM]{lavaan2RAM}}.}
\item{data}{A data frame or matrix of data. Alternatively, users may
  provide a \code{\link[OpenMx]{mxData}} object.}
\item{intervals.type}{Either \code{z} (default if missing) or
  \code{LB}. If it is \code{z}, it calculates the 95\% confidence
  intervals (CIs) based on the estimated standard error. If it
  is \code{LB}, it calculates the 95\% likelihood-based CIs on the
  parameter estimates.}
\item{var.names}{If it is \code{NULL} (the default), the
  \code{var.names} are extracted from \code{Fmatrix$values}.}
\item{mxModel.Args}{A list of arguments passed to \code{\link[OpenMx]{mxModel}}.}
\item{run}{Logical. If \code{FALSE}, only return the mx model without
  running the analysis.}
\item{mxTryHard}{If \code{TRUE}, call \code{mxTryHard} to run the
  analysis. If \code{FALSE}, call \code{mxRun} to run the analysis.}
\item{silent}{Logical. An argument is passed to \code{\link[OpenMx]{mxRun}}}

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

## A regression model
model <- "y ~ x   # Regress y on x
          y ~ 1   # Intercept of y
          x ~ 1   # Mean of x"
plot(model)

RAM <- lavaan2RAM(model, obs.variables=c("y", "x"))

my.fit <- create.mxModel(RAM=RAM, data=my.df)
summary(my.fit)

## A meta-analysis
model <- "f =~ 1*yi
          f ~ mu*1          ## Average effect
          f ~~ tau2*f       ## Heterogeneity variance
          yi ~~ data.vi*yi  ## Known sampling variance"
plot(model)

## Do not standardize the latent variable (f): std.lv=FALSE 
RAM <- lavaan2RAM(model, obs.variables="yi", std.lv=FALSE)

## Use likelihood-based CI 
my.fit <- create.mxModel(RAM=RAM, data=Hox02, intervals="LB")
summary(my.fit)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
