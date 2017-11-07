\name{uniR2}
\alias{uniR2mx}
\alias{uniR2lavaan}
\title{Second Stage analysis of the univariate R (uniR) approach
}
\description{It conducts the second stage analysis of the uniR analysis
  by fitting structural equation models on the average correlation matrix. 
}
\usage{
uniR2mx(x, Amatrix = NULL, Smatrix = NULL, Fmatrix = NULL,
        model.name=NULL, suppressWarnings=TRUE, silent=TRUE,
        run=TRUE, ...)
uniR2lavaan(x, model, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{uniR1} from
	\code{\link[metaSEM]{uniR1}}.
}
  \item{Amatrix}{An asymmetric matrix in the RAM specification with
	\code{\link[OpenMx]{MxMatrix-class}}. If it is a matrix, it will be
	converted into \code{\link[OpenMx]{MxMatrix-class}} by the \code{as.mxMatrix} function.
}
  \item{Smatrix}{A symmetric matrix in the RAM specification with
	\code{\link[OpenMx]{MxMatrix-class}}. If it is a matrix, it will be
	converted into \code{\link[OpenMx]{MxMatrix-class}} by the \code{as.mxMatrix} function.
}
  \item{Fmatrix}{A filter matrix in the RAM specification with
	\code{\link[OpenMx]{MxMatrix-class}}. If it is \code{NULL} (the default), an
	identity matrix with the same dimensions of \code{Cov} will be
	created. If it is a matrix, it will be converted into
	\code{\link[OpenMx]{MxMatrix-class}} by the \code{as.mxMatrix} function. It is not required when there is no latent variable.
}
  \item{model.name}{A string for the model name in
	\code{\link[OpenMx]{mxModel}}. If it is missing, the default is
	"UniR2".
}
  \item{suppressWarnings}{Logical. If \code{TRUE}, warnings are
	suppressed. Argument to be passed to \code{\link[OpenMx]{mxRun}}.
}
  \item{silent}{Logical. Argument to be passed to \code{\link[OpenMx]{mxRun}}}
  \item{run}{Logical. If \code{FALSE}, only return the mx model without
	running the analysis.} 
  \item{model}{A model specified using lavaan syntax see \code{\link[lavaan]{model.syntax}}}  
  \item{\dots}{Further arguments to be passed to either
	\code{\link[OpenMx]{mxRun}} or \code{\link[lavaan]{sem}. For
	  \code{\link[lavaan]{sem}}, \code{fixed.x=FALSE} is passed automatically.}
  }
}

\value{A fitted object returned from \code{\link[OpenMx]{mxRun}} or \code{\link[lavaan]{sem}}.
}
\details{
  This function implements the univariate r approach proposed by Viswesvaran
  and Ones (1995) to conduct meta-analytic structural equation modeling
  (MASEM). It treats the average correlation matrix as if it was a
  covariance matrix in fitting structural equation models. The harmonic
  mean of the sample sizes in combining correlation coefficients is used
  as the sample size in fitting structural equation models. It is
  included in this package for research interests. The two-stage structural equation
  modeling (TSSEM) approach is preferred (e.g., Cheung, 2015; Cheung & Chan, 2005).
}
\references{
  Cheung, M. W.-L. (2015). \emph{Meta-analysis: A structural equation
	modeling approach}. Chichester, West Sussex: John Wiley & Sons, Inc.

  Cheung, M. W.-L., & Chan, W. (2005). Meta-analytic structural equation
  modeling: A two-stage approach. \emph{Psychological Methods},
  \bold{10}, 40-64.

  Viswesvaran, C., & Ones, D. S. (1995). Theory testing: Combining
  psychometric meta-analysis and structural equations
  modeling. \emph{Personnel Psychology}, \bold{48}, 865-885. 
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\seealso{ \code{\link[metaSEM]{uniR1}}, \code{\link[metaSEM]{lavaan2RAM}},
    \code{\link[metaSEM]{Becker09}}
}
\keyword{uniR}
