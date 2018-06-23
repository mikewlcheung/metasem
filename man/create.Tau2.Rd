\name{create.Tau2}
\alias{create.Tau2}
\title{Create variance component of the heterogeneity of the random effects
}
\description{It creates variance component of the heterogeneity of the
  random effects by decomposing the variance component into matrices of
  correlation and standard deviations.  
}
\usage{
create.Tau2(RAM, no.var, RE.type = c("Diag", "Symm", "Zero"),
            Transform = c("expLog", "sqSD"), RE.startvalues=0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{RAM}{The RAM model for testing. \code{no.var} is calculated from it.}  
\item{no.var}{If \code{RAM} is missing, the user must have to specify
  the \code{no.var} argument. It represents the number of variables (random effects).}
\item{RE.type}{Either \code{"Diag"}, \code{"Symm"} or \code{"Zero"}. If
  it is\code{"Diag"} (the default if missing), a diagonal matrix is used
  for the random effects meaning that the random effects are
  independent. If it is \code{"Symm"}, a symmetric matrix is used for the random effects on the covariances
	among the correlation (or covariance) vectors. If it is
	\code{"Zero"}, a zero matrix is assumed on the variance component of
  the random effects.}
\item{Transform}{Either \code{"expLog"} or \code{"sqSD"}. If it is
  \code{"expLog"}, the variances are estimated by applying a log and exp
  transformation. If it is \code{"sqSD"}, the variances are estimated by
  applying a square on the SD. The transformation may improve the
  estimation when the heterogeneity is small or close to zero.}   
\item{RE.startvalues}{Starting values for the variances.}
}

\value{A list of \code{MxMatrix-class}. The variance component is
  computed in \code{Tau2}.}

\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\seealso{ \code{\link[metaSEM]{osmasem}},
  \code{\link[metaSEM]{create.V}}, \code{\link[metaSEM]{create.vechsR}}
}

\examples{
\dontrun{
T0 <- create.Tau2(no.var=6, RE.type="Diag", Transform="expLog", RE.startvalues=0.05)
T0    
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ osmasem }
