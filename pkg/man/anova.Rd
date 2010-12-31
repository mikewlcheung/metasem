\name{anova}
\alias{anova.wls}
\alias{anova.meta}
\alias{anova.reml}
\title{Compare Models with Likelihood Ratio Statitics
}
\description{It compares models with likelihood ratio statistic from either \code{wls}, \code{meta} or
  \code{reml} objects via \code{\link[OpenMx]{mxCompare}}.
}
\usage{
\method{anova}{wls}(object, \dots, all=FALSE)
\method{anova}{meta}(object, \dots, all=FALSE)
\method{anova}{reml}(object, \dots, all=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object or a list of objects returned from either class
	\code{\link[metaSEM]{wls}}, class \code{\link[metaSEM]{meta}} or
	class \code{\link[metaSEM]{reml}}. It will be passed to the
	\code{base} argument in \code{\link[OpenMx]{mxCompare}}.
  }
  \item{\dots}{An object or a list of objects returned from either class
	\code{\link[metaSEM]{wls}}, class \code{\link[metaSEM]{meta}} or
	class \code{\link[metaSEM]{reml}}. It will be passed to the
	\code{comparison} argument in \code{\link[OpenMx]{mxCompare}}.
  }  
  \item{all}{A boolean value on whether to compare all bases with all
	comparisons. It will be passed to the \code{all} argument in
	\code{\link[OpenMx]{mxCompare}}.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{A table of comparisons between the models in base and comparison.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\note{Special care should be taken to make sure that the models between
  \code{base} and \code{comparison} are nested. One common mistake is
  comparing a model without predictor and a model with predictor. Since
  the predictors are also involved in the estimation in \code{metaSEM},
  these two models are not nested. The correct way to compare them is to
  fix the regression coefficient of one model at zero while the
  coefficient in the other model is free.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## Test the significance of a predictor with likelihood ratio test
## Model0: Regression coefficient is fixed at 0
model0 <- with( Hox02, meta(y=yi, v=vi, x=weeks, coeff.constraints=matrix(0),
                            model.name="No predictor") )

## Model1: Regression coefficient is free
model1 <- with( Hox02, meta(y=yi, v=vi, x=weeks,
                            model.name="One predictor") )
anova(model1, model0) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

