\name{lavaan2RAM}
\alias{lavaan2RAM}
\title{
Convert \code{lavaan} models to RAM models
}
\description{
It converts models specified in \code{lavaan} model syntax to RAM models.
}
\usage{
lavaan2RAM(model, obs.variables = NULL, A.notation = "ON",
           S.notation = "WITH", M.notation = "mean", auto.var = TRUE,
           std.lv = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{A character string of model using the lavaan model syntax.
}
\item{obs.variables}{A character vector of the observed variables. The
  observed variables in the RAM specification will follow the order
  specified in \code{obs.variables}. It is important to check whether the
  order of the observed variables matches the order in the dataset. 
}
\item{A.notation}{A character string to be used in the A matrix if the labels
  are not included in the lavaan model. For example, the label will be "yONx" for
  regressing "y" on "x".
}
\item{S.notation}{A character string to be used in the S matrix if the labels
  are not included in the lavaan model. For example, the label will be "yWITHx" for
  the covariance between "y" with "x" and "yWITHy" for the (error) variance
  of "y". 
}
\item{M.notation}{A character string to be used in the M matrix if the labels
  are not included in the lavaan model. For example, the label will be "ymean" for
  the mean of "y" if \code{M.notation="mean"}. 
}
\item{auto.var}{Logical. If \code{TRUE}, the residual variances and the
  variances of exogenous latent variables are included in the model and
  set free. See \code{\link[lavaan]{lavaanify}}.
}
\item{std.lv}{Logical. If \code{TRUE}, the metric of each latent
  variable is determined by fixing their variances to 1.0. If FALSE, the
  metric of each latent variable is determined by fixing the factor
  loading of the first indicator to 1.0. See \code{\link[lavaan]{lavaanify}}.
}
\item{\dots}{Further arguments to be passed to \code{\link[lavaan]{lavaanify}}}. Please note that \code{fixed.x} is set at \code{FALSE}. Thus it cannot be passed to \dots.
}
\details{
It uses the \code{\link[lavaan]{lavaanify}} to do the conversion.
}
\value{
  A list of RAM specification with \code{A}, \code{S}, \code{F}, and \code{M} matrices.
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}

\seealso{\code{\link[semPlot]{ramModel}},
  \code{\link[metaSEM]{Becker92}}, \code{\link[metaSEM]{Becker09}},
  \code{\link[metaSEM]{Digman97}}, \code{\link[metaSEM]{Hunter83}},
  \code{\link[metaSEM]{as.mxMatrix}},
  \code{\link[metaSEM]{checkRAM}}
}
\examples{
## Regression model on correlation matrix
model1 <- "## y is modelled by x1, x2, and x3
           y ~ x1 + x2 + x3
           ## Fix the independent variables at 1
           x1 ~~ 1*x1
           x2 ~~ 1*x2
           x3 ~~ 1*x3
           ## Declare the correlations among the independent variables
           x1 ~~ x2
           x1 ~~ x3
           x2 ~~ x3"

## Compare the arrangements of variables with and without
## specifying the obs.variables arguments. 
lavaan2RAM(model1, obs.variables=c("y", "x1", "x2", "x3"))

## Two-factor CFA model
model2 <- "f1 =~ x1 + x2 + x3
           f2 =~ x4 + x5 + x6
           ## Declare the correlation between f1 and f2
           ## and label it with cor_f1f2
           f1 ~~ cor_f1f2*f2"

lavaan2RAM(model2)

## Regression model with the mean structure
model3 <- "y ~ x
           ## Intercept of y
           y ~ 1
           ## Mean of x
           x ~ 1"

lavaan2RAM(model3)
}
\keyword{ methods }
\keyword{ tssem }% __ONLY ONE__ keyword per line
