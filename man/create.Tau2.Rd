\name{create.Tau2}
\alias{create.Tau2}
\title{Create a variance component of the heterogeneity of the random effects
}
\description{It creates variance component of the heterogeneity of the
  random effects by decomposing the variance component into matrices of
  correlation and standard deviations.  
}
\usage{
create.Tau2(RAM, no.var, Tau1.labels=seq(no.var),
            RE.type = c("Diag", "Symm", "Zero", "User"),
            RE.User=NULL, Transform = c("expLog", "sqSD"),
            RE.startvalues=0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{RAM}{The RAM model for testing. \code{no.var} is calculated from it.}  
\item{no.var}{If \code{RAM} is missing, the user has to specify
  the \code{no.var} argument. It represents the \code{no.var} by \code{no.var} of the random effects).}
\item{Tau1.labels}{Parameter labels in \code{Tau1}. The default is
  \code{Tau1_1}, \code{Tau1_2}, etc.
}
\item{RE.type}{Either \code{"Diag"}, \code{"Symm"}, \code{"Zero"} or \code{"User"}. If
  it is\code{"Diag"} (the default if missing), a diagonal matrix is used
  for the random effects meaning that the random effects are
  independent. If it is \code{"Symm"}, a symmetric matrix is used for the random effects on the covariances
	among the correlation (or covariance) vectors. If it is
	\code{"Zero"}, a zero matrix is assumed on the variance component of
  the random effects. If it is \code{"User"}, users have to specify the
  \code{RE.User} argument.}
\item{RE.User}{It represents the \code{no.var} by \code{no.var} symmetric matrix
  of \code{TRUE} or \code{FALSE} for the variance component. If the
  elements are \code{FALSE}, they are fixed at 0.}
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
T0 <- create.Tau2(no.var=4, RE.type="Diag", Transform="expLog", RE.startvalues=0.05)
T0

T1 <- create.Tau2(no.var=4, Tau1.labels=c("a", "b", "c", "d"))
T1
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ osmasem }
