\name{deltaMethod}
\alias{deltaMethod}
\title{Approximate Variance or Covariance Matrix of Functions Using the Delta Method
}
\description{It calculates the approximate variance or covariance matrix
  of functions using the delta method. 
}
\usage{
deltaMethod(fn, variables, Mean, Cov)
}
\arguments{
  \item{fn}{A list of functions in characters.  
}  
  \item{variables}{Variables used in \code{fn}, \code{Mean}, and \code{Cov}
}
  \item{Mean}{The sample mean or the parameter estimates
}
  \item{Cov}{The sample covariance matrix or the sampling covariance of
	the paramter estimates
}	
}
\note{This function is modified from the \code{deltamethod()} in the
  \code{msm} package}
\value{The approximate variance or the covariance matrix of \code{fn}
}

\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
\dontrun{
## Example taken from http://lavaan.ugent.be/tutorial/mediation.html
## Generate artificial data to estimate the indirect effect
set.seed(1234)
X <- rnorm(100)
M <- 0.5*X + rnorm(100)
Y <- 0.7*M + rnorm(100)
Data <- data.frame(X=X, Y=Y, M=M)

model <- '# direct effect
          Y ~ c*X
          # mediator
          M ~ a*X
          Y ~ b*M
          # indirect effect (a*b)
          ab := a*b
          # total effect
          total := c + (a*b)'
fit <- sem(model, data = Data)
summary(fit)
}


## Approximate standard errors of the indirect and total effects using the delta method.
## The first and second functions are the direct and indirect effects, respectively.
fn <- list("a*b", "a*b+c")
variables <- c("a", "b", "c")

## Results from the above analysis
Mean <- c(0.47391503, 0.78831930, 0.03635467)

Cov <- structure(c(0.0105540501465901, -5.36944694410567e-19,
4.78727359106282e-18, 3.56339238409499e-19, 0.00852258573549811,
-0.00403898148221173, 1.60078511618195e-18, -0.00403898148221173,
0.0109089137570649), .Dim = c(3L, 3L), .Dimnames = list(c("a", "b", "c"), c("a", 
"b", "c")))

ind.Cov <- deltaMethod(fn=fn, variables=variables, Mean=Mean, Cov=Cov)

## Standard errors
sqrt(diag(ind.Cov))


## Approximate sampling variance of log(p/(1-p))
## Number of success
x <- 20
## Total number of trials
n <- 100
## Estimated proportion of success
p <- x/n
## Sampling variance of p
v <- p*(1-p)/n

## logit(p)
fn <- "log(p/(1-p))"
variables <- c("p")

## Sampling variance of logit(p)
deltaMethod(fn=fn, variables=variables, Mean=p, Cov=v)
}
\keyword{utilities}
