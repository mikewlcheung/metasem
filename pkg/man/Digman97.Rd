\name{Digman97}
\alias{Digman97}
\docType{data}
\title{Factor Correlation Matrices of Big Five Model from Digman (1997)}
\description{
The data set includes fourteen studies of the factor correlation matrices of the Five-Factor Model of personality reported by Digman (1997). 
}
\usage{data(Digman97)}
\format{
}
\details{
A list of data with the following structure:
\describe{
  \item{data}{A list of 14 studies of correlation matrices. The
	variables are \emph{Extraversion} (E), \emph{Agreeableness} (A), \emph{Conscientiousness}
	(C), \emph{Emotional Stability} (ES) and \emph{Intellect} (I)}
  \item{n}{A vector of sample sizes}
  \item{cluster}{Types of participants of the studies}
}
}
\source{
Digman, J.M. (1997). Higher-order factors of the Big Five. \emph{Journal of Personality and Social Psychology}, \bold{73}, 1246-1256.
}
\references{
 Cheung, M. W. L., & Chan, W. (2005). Classifying correlation matrices into relatively homogeneous subgroups: A cluster analytic approach. \emph{Educational and Psychological Measurement}, \bold{65}, 954-979.
}
\examples{
Digman97

##### Example of Fixed-effects TSSEM
fixed1 <- tssem1(Digman97$data, Digman97$n, method="FEM")
summary(fixed1)

## Factor correlation matrix
Phi <- mxMatrix("Stand", ncol=2, nrow=2, value=.2, free=TRUE, name="Phi")
## Factor loadings
Lambda <- mxMatrix("Full", ncol=2, nrow=5, value=c(0,.3,.3,.3,0,.3,0,0,0,.3),
free=c(FALSE,TRUE,TRUE,TRUE,FALSE,TRUE,FALSE,FALSE,FALSE,TRUE),
name="Lambda")

## Model to be fitted
impliedR <- mxAlgebra(Lambda \%&\% Phi, name="impliedR")

fixed2 <- tssem2(fixed1, impliedS=impliedR, matrices=c(Phi, Lambda),
                  model.name="TSSEM2 Digman97")
summary(fixed2)

##### Example of Fixed-effects TSSEM with several clusters
fixed1.cluster <- tssem1(Digman97$data, Digman97$n, cluster=Digman97$cluster)
summary(fixed1.cluster)

fixed2.cluster <- tssem2(fixed1.cluster, impliedS=impliedR, matrices=c(Phi, Lambda))
summary(fixed2.cluster)

##### Example of Random-effects TSSEM with random effects on the diagonals
random1 <- tssem1(Digman97$data, Digman97$n, method="REM", RE.diag.only=TRUE)
summary(random1)

random2 <- tssem2(random1, impliedS=impliedR, matrices=c(Phi, Lambda))
summary(random2)
}
\keyword{datasets}
