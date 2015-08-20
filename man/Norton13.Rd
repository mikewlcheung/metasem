\name{Norton13}
\alias{Norton13}
\docType{data}
\title{Studies on the Hospital Anxiety and Depression Scale Reported by Norton et al. (2013)
}
\description{
  The data set includes 28 studies on 14 items measuring the
  Hospital Anxiety and Depression Scale (HADS) Reported by Norton et
  al. (2013).
}
\usage{data(Norton13)}

\details{
      The variables are:
  \describe{
  \item{data}{A list of 28 studies of correlation matrices. The
      variables are 14 items (x1 to x14) measuring HADS.}
  \item{n}{A vector of sample sizes}
  \item{population}{A vector of population of the data}
  \item{group}{A vector of classification into \emph{patients}
      vs. \emph{non-patients} based on population}
  }
}
\source{
Norton, S., Cosco, T., Doyle, F., Done, J., & Sacker, A. (2013). The Hospital Anxiety and Depression Scale: A meta confirmatory factor analysis. \emph{Journal of Psychosomatic Research}, \emph{74}(1), 74-81.
}
\references{
  Jak, S., & Cheung, M. W.-L. (2015). Addressing heterogeneity in meta-analytic structural equation modeling using subgroup
analysis. \emph{Manuscript in preparation}.
}
\examples{
\dontrun{
## Fixed-effects model: First stage analysis
fixed1 <- tssem1(my.df = Norton13$data, n = Norton13$n, method = "FEM")
summary(fixed1)

## Names of the observed variables
varnames <- paste("x", 1:14, sep="")

## Names of the latent variables
factornames <- c("Dist","Anx","Dep")

## Vector to indicate whether variable is observed (1) or latent (0)
F <- create.Fmatrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0), name="F")

## Matrix with factor loadings, later included in matrix A
lambda <-matrix(
		c("0.1*L1_1","0.1*L1_2","0",
		  "0.1*L2_1","0","0.1*L2_3",
		  "0.1*L3_1","0.1*L3_2","0",
		  "0.1*L4_1","0","0.1*L4_3",
		  "0.1*L5_1","0.1*L5_2","0",
		  "0.1*L6_1","0","0.1*L6_3",
		  "0.1*L7_1","0.1*L7_2","0",
		  "0.1*L8_1","0","0.1*L8_3",
		  "0.1*L9_1","0.1*L9_2","0",
		  "0.1*L10_1","0","0.1*L10_3",
		  "0.1*L11_1","0.1*L11_2","0",
		  "0.1*L12_1","0","0.1*L12_3",
		  "0.1*L13_1","0.1*L13_2","0",
		  "0.1*L14_1","0","0.1*L14_3"), 
	   	nrow=14, ncol=3, byrow = TRUE)

A <- rbind(cbind(matrix(0,ncol=14,nrow=14),lambda),
           matrix(0, nrow=3, ncol=17))

## This step is not necessary but it is useful for inspecting the model.
dimnames(A) <- list(c(varnames,factornames), c(varnames,factornames))

## Display A
A

## Matrix with residual variances
psi <- matrix(0, nrow = 14, ncol = 14)
diag(psi) <- paste("0.1*","e",1:14,sep="")

# Matrix with factor variances and covariances
phi <- matrix(c(1,0,0,
                0,1,0,
                0,0,1), nrow = 3, ncol = 3)

# Psi and phi are combined in matrix S
S <- bdiagMat(list(psi, phi))

## This step is not necessary but it is useful for inspecting the model.
dimnames(S) <- list(c(varnames,factornames),c(varnames,factornames))

## Display S
S

## Fixed-effects model: Second stage analysis
fixed2 <- tssem2(fixed1, Amatrix=A, Smatrix=S, diag.constraints=TRUE, 
                 intervals="LB")
summary(fixed2)

## Load the library
library("semPlot")

## Convert the model to semPlotModel object
my.plot <- meta2semPlot(fixed2)

## Plot the parameter estimates
semPaths(my.plot, whatLabels="est", nCharNodes=10, color="green")
}
}
\keyword{datasets}
