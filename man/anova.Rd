\name{anova}
\alias{anova.wls}
\alias{anova.meta}
\alias{anova.meta3LFIML}
\alias{anova.reml}
\alias{anova.osmasem}
\title{Compare Nested Models with Likelihood Ratio Statistic
}
\description{It compares nested models with the likelihood ratio
  statistic from various objects. It is a wrapper of \code{\link[OpenMx]{mxCompare}}.
}
\usage{
\method{anova}{wls}(object, \dots, all=FALSE)
\method{anova}{meta}(object, \dots, all=FALSE)
\method{anova}{meta3LFIML}(object, \dots, all=FALSE)
\method{anova}{reml}(object, \dots, all=FALSE)
\method{anova}{osmasem}(object, \dots, all=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object or a list of objects of various classes. It will be passed to the
	\code{base} argument in \code{\link[OpenMx]{mxCompare}}.
  }
  \item{\dots}{An object or a list of objects of various classes. It will be passed to the
	\code{comparison} argument in \code{\link[OpenMx]{mxCompare}}.
  }  
  \item{all}{A Boolean value on whether to compare all bases with all
	comparisons. It will be passed to the \code{all} argument in
	\code{\link[OpenMx]{mxCompare}}.}
}
% \details{Special care has to be taken to make sure that the models being
%   compared in \code{base} and \code{comparison} are nested. One common mistake is
%   comparing a model without predictor and a model with predictors in \code{\link[metaSEM]{meta}}. Since
%   the parameters in the predictors, e.g., means and variances, are also estimated in \code{\link[metaSEM]{meta}},
%   these two models are not nested. The correct way to compare them is to
%   fix the regression coefficients of one model at zero while the
%   coefficients in the other model are free (see the example). If only one parameter is tested, an alternative (and easier) approach
%   is to request the likelihood-based CI directly with the
%   \code{intervals.type="LB"} argument in \code{\link[metaSEM]{meta}}.
%
% }
\value{A table of comparisons between the models in base and comparison.
}

\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\note{When the objects are class \code{\link[metaSEM]{wls}}, the degrees
  of freedom in the base and comparison models are incorrect, while the degrees of
  freedom of the difference between them is correct. If users want to
  obtain the correct degrees of freedom in the base and comparison
  models, they may individually apply the \code{\link[metaSEM]{summary}} function on
  the base and comparison models.
}
  
\examples{
## Test the significance of a predictor with likelihood ratio test
## Model0: No predictor
model0 <- meta(y=yi, v=vi, data=Hox02, model.name="No predictor")

## Model1: With a predictor
model1 <- meta(y=yi, v=vi, x=weeks, data=Hox02, model.name="One predictor")

## Compare these two models
anova(model1, model0) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

