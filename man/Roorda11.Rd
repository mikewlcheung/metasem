\name{Roorda11}
\alias{Roorda11}
\docType{data}
\title{Studies on Students' School Engagement and Achievement Reported
  by Roorda et al. (2011)
}
\description{
The data set includes 45 studies on the influence of affective teacher–student relationships on students’ school engagement and
achievement reported by Roorda et al. (2011).
}
\usage{data(Roorda11)}

\details{
      The variables are:
  \describe{
  \item{data}{A list of 45 studies of correlation matrices. The
      variables are \emph{pos} (positive teacher-student relations),
      \emph{neg} (negative teacher-student relations), \emph{enga}
      (student engagement), and \emph{achiev} (student achievement).}
  \item{n}{A vector of sample sizes}
  \item{SES}{A vector of average socio-economic status (SES) of the samples}
  }
}
\source{
  Roorda, D. L., Koomen, H. M. Y., Spilt, J. L., & Oort, F. J. (2011). The influence of affective teacher–student relationships on students' school engagement and achievement a meta-analytic approach. \emph{Review of Educational Research}, \emph{81}(4), 493-529.
}
\references{
  Jak, S., & Cheung, M. W.-L. (2015). Addressing heterogeneity in meta-analytic structural equation modeling using subgroup
analysis. \emph{Manuscript in preparation}.
}
\examples{
\dontrun{

## Random-effects model: First stage analysis
random1 <- tssem1(my.df = Roorda11$data, n = Roorda11$n, method = "REM", RE.type = "Diag")
summary(random1)

varnames <- c("pos", "neg", "enga", "achiev")

## Prepare a regression model using create.mxMatrix()
A <- create.mxMatrix(c(0,0,0,0,
                       0,0,0,0,
                       "0.1*b31","0.1*b32",0,0,
                       0,0,"0.1*b43",0),
                     type = "Full", nrow = 4, ncol = 4, byrow = TRUE,
                     name = "A", as.mxMatrix = FALSE)

## This step is not necessary but it is useful for inspecting the model.
dimnames(A) <- list(varnames, varnames)
A

S <- create.mxMatrix(c(1,
                       ".5*p21",1,
                       0,0,"0.6*p33",
                       0,0,0,"0.6*p44"), 
                     type="Symm", byrow = TRUE,
                     name="S", as.mxMatrix = FALSE)

## This step is not necessary but it is useful for inspecting the model.
dimnames(S) <- list(varnames, varnames)
S

## Random-effects model: Second stage analysis
random2 <- tssem2(random1, Amatrix=A, Smatrix=S, diag.constraints=TRUE, 
                  intervals="LB")
summary(random2)

## Load the library
library("semPlot")

## Convert the model to semPlotModel object
my.plot <- meta2semPlot(random2)

## Plot the parameter estimates
semPaths(my.plot, whatLabels="est", nCharNodes=10, color="green")
}
}
\keyword{datasets}
