\name{issp05}
\alias{issp05}
\docType{data}
\title{Data Set from ISSP (2005) 
}
\description{Thirty-two covariance matrices on work-related attitudes
  were extracted from the International Social Survey Programme 2005:
  Work Orientation III (ISSP, 2005). Seven variables were selected for
  demonstration purposes. They were grouped into three constructs:
  \emph{Importance of Job Prospects} measured by job security (JP1),
  high income (JP2), and opportunity for advancement (JP3); \emph{Importance
  of Job Autonomy} measured by work independently (JA1) and decide time
  of work (JA2); and \emph{Importance of Contributions to Society} measured by help other
  people (CS1) and a job useful to society (CS2).
}
\usage{data(issp05)}
\format{
}
\details{
A list of data with the following structure:
\describe{
  \item{data}{A list of 32 covariance matrices }
  \item{n}{A vector of sample sizes}
}
}
\source{ISSP Research Group (2007): \emph{International Social Survey Programme 2005: Work Orientation III (ISSP 2005)}. GESIS Data Archive, Cologne. ZA4350 Data file Version 1.0.0, doi:10.4232/1.4350}

\references{
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\seealso{ \code{\link[metaSEM]{issp89}}
}

\examples{
issp05

#### Fixed-effects TSSEM
fixed1 <- tssem1(issp05$data, issp05$n, method="FEM")
summary(fixed1)

Phi <- mxMatrix("Stand", ncol=3, nrow=3, value=.2, free=TRUE, name="Phi")
Lambda <- mxMatrix("Full", ncol=3, nrow=7,
                   value=c(rep(0.2,3), rep(0,7), rep(0.2,2), rep(0,7), rep(0.2,2)),
                   free=c(rep(TRUE,3), rep(FALSE,7), rep(TRUE,2), rep(FALSE,7), rep(TRUE,2)),
                   name="Lambda")

## Model to be fitted
impliedR <- mxAlgebra(Lambda \%&\% Phi, name="impliedR")

fixed2 <- tssem2(fixed1, impliedS=impliedR, matrices=c(Phi, Lambda),
                 model.name="TSSEM2 ISSP05")
summary(fixed2)
}
\keyword{datasets}
