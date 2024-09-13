\name{sem}
\alias{sem}
\alias{create.mxModel}
\title{Fit a structural equation model using OpenMx
}
\description{It fits a structural equation model by creating a mxModel from a RAM object.  
}
\usage{
## Depreciated in the future
create.mxModel(model.name="sem", RAM=NULL, data=NULL,
               Cov=NULL, means=NULL, numObs,
               intervals.type=c("z", "LB"), startvalues=NULL,
               replace.constraints=FALSE, mxModel.Args=NULL,
               run=TRUE, silent=TRUE, ...)
sem(model.name="sem", RAM=NULL, data=NULL,
    Cov=NULL, means=NULL, numObs,
    intervals.type=c("z", "LB"), startvalues=NULL,
    lbound=NULL, ubound=NULL, replace.constraints=FALSE,
    mxModel.Args=NULL, run=TRUE, silent=TRUE, ...)
}
\arguments{
  \item{model.name}{A string for the model name in \code{\link[OpenMx]{mxModel}}. 
}
\item{RAM}{A RAM object including a list of matrices of the model
  returned from \code{\link[metaSEM]{lavaan2RAM}}.}
\item{data}{A data frame or matrix of data.}
\item{Cov}{A covariance matrix may also be used if \code{data}==NULL.}
\item{means}{A named vector of means (options) if \code{Cov} is used.}
\item{numObs}{If \code{Cov} is used, a sample size must be provided.}
\item{intervals.type}{Either \code{z} (default if missing) or
  \code{LB}. If it is \code{z}, it calculates the 95\% confidence
  intervals (CIs) based on the estimated standard error. If it
  is \code{LB}, it calculates the 95\% likelihood-based CIs on the
  parameter estimates.}
\item{startvalues}{A list of named starting values of the free parameters, e.g., list(a=1, b=2)}
\item{lbound}{A list of lower bound of the free parameters. If it is
  not provided, all free parameters are assumed \code{NA}.}
\item{ubound}{A list of upper bound of the free parameters. If it is
  not provided, all free parameters are assumed \code{NA}.}  
\item{replace.constraints}{Logical. If \code{TRUE}, the parameters on the left hand side will be replaced by the constraints on the right hand side. That is, the parameters on the left hand side are no longer parameters in the model.}
\item{mxModel.Args}{A list of arguments passed to \code{\link[OpenMx]{mxModel}}.}
\item{run}{Logical. If \code{FALSE}, only return the mx model without
  running the analysis.}
\item{silent}{Logical. An argument is passed to either
  \code{\link[OpenMx]{mxRun}} or \code{\link[OpenMx]{mxTryHard}}}
\item{\dots}{Further arguments will be passed to either
  \code{\link[OpenMx]{mxRun}} or \code{\link[OpenMx]{mxTryHard}}}  
}
\value{An object of class \code{mxsem}}
\note{when there are constraints with \code{replace.constraints=TRUE} and \code{intervals.type="LB"}, it returns an error because some parameters in the model are replaced with the new parameters in the constraints. However, the names of these new parameters are not captured in the CI object. Users are advised to use \code{intervals.type="z"} before it is fixed.
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}

\examples{
\donttest{
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

my.fit <- sem(RAM=RAM, data=my.df)
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
my.fit <- sem(RAM=RAM, data=Hox02, intervals="LB")
summary(my.fit)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
