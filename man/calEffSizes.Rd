\name{calEffSizes}
\alias{calEffSizes}
\title{Calculate Effect Sizes using lavaan Models
}\description{It calculates effect sizes with Delta Method by formulating the effect sizes
  as functions of SEM in lavaan.
}
\usage{
calEffSizes(model, data=NULL, n, Cov, Mean=NULL, group=NULL, lavaan.output=FALSE, ...)
}
\arguments{
  \item{model}{A lavaan model. Effect sizes are defined as functions of
	SEM parameters with \code{:=}.}
  \item{data}{A data frame of the observed variables. If it is
	\code{NULL}, summary statistics are required.}
  \item{n}{Sample sizes}
  \item{Cov}{A covariance matrix or a list of covariance matrices.}
  \item{Mean}{Optional sample means.}
  \item{group}{A character of the variable name in the data frame defining the groups in a multiple group analysis.}
  \item{lavaan.output}{If \code{TRUE}, it returns the fitted object
	instead of the effect sizes and their sampling covariance matrix.}
  \item{\dots}{Further arguments passed to \code{\link[lavaan]{sem}}.}
}

\value{Effect sizes and their sampling covariance matrix or a lavaan
  fitted object.
}
\note{The input matrices are treated as covariance matrices unless there
  are explicit constraints in the model.
}
\references{
  Cheung, M. W.-L. (2015). \emph{Meta-analysis: A structural equation
	modeling approach}. Chichester, West Sussex: John Wiley & Sons, Inc.

  Cheung, M. W.-L. (2018). Computing multivariate effect sizes and their sampling covariance matrices with structural equation modeling: Theory, examples, and computer simulations. \emph{Frontiers in Psychology}, \bold{9}(1387). https://doi.org/10.3389/fpsyg.2018.01387
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link[metaSEM]{smdMES}}, \code{\link[metaSEM]{smdMTS}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\dontrun{ 
## Select ATT, Bi, and BEH 
obs.vars <- c("BEH", "BI", "ATT")

## Select one study from Cooke16 for illustration
my.cor <- Cooke16$data[[4]][obs.vars, obs.vars]
my.n  <- Cooke16$n[4]

## Effect sizes: indirect effect and direct effect
model <- "BEH ~ c*ATT + b*BI
          BI ~ a*ATT
          ## Indirect effect
          Ind := a*b
          Dir := c"

calEffSizes(model=model, n=my.n, Cov=my.cor, lavaan.output=FALSE)

## Return the lavaan fitted model
fit <- calEffSizes(model=model, n=my.n, Cov=my.cor, lavaan.output=TRUE)
lavaan::summary(fit)
    
lavaan::parameterestimates(fit)  
}
}
\keyword{meta-analysis }
