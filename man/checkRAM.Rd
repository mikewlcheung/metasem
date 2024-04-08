\name{checkRAM}
\alias{checkRAM}
\title{Check the correctness of the RAM formulation
}
\description{It provides simple checks on the correctness of the RAM formulation.
}
\usage{
checkRAM(Amatrix, Smatrix, cor.analysis=TRUE)
}
\arguments{
 \item{Amatrix}{An asymmetric matrix in the RAM specification with
   \code{\link[OpenMx]{MxMatrix-class}}. If it is a matrix, it will be converted into \code{\link[OpenMx]{MxMatrix-class}} by the \code{as.mxMatrix} function.
}
  \item{Smatrix}{A symmetric matrix in the RAM specification with
	\code{\link[OpenMx]{MxMatrix-class}}. If it is a matrix, it will be converted into \code{\link[OpenMx]{MxMatrix-class}} by the \code{as.mxMatrix} function.
  }
  \item{cor.analysis}{Logical. Analysis of correlation or covariance
	structure. There are additional checks for cor.analysis=\code{TRUE}.}
}
\value{It returns silently if no error has been detected; otherwise, it
  returns a warning message.
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}

\seealso{ \code{\link[metaSEM]{as.mxMatrix}},
  \code{\link[metaSEM]{lavaan2RAM}}
}
\examples{
\donttest{ 
## Digman97 example
model1 <- "## Factor loadings
           Alpha=~A+C+ES
           Beta=~E+I
           ## Factor correlation
           Alpha~~Beta"

RAM1 <- lavaan2RAM(model1, obs.variables=c("A","C","ES","E","I"),
                   A.notation="on", S.notation="with")
RAM1

## The model is okay.    
checkRAM(Amatrix=RAM1$A, Smatrix=RAM1$S)

## Hunter83 example    
model2 <- "## Regression paths
           Job_knowledge ~ A2J*Ability
           Work_sample ~ A2W*Ability + J2W*Job_knowledge
           Supervisor ~ J2S*Job_knowledge + W2S*Work_sample

           ## Fix the variance of Ability at 1
           Ability ~~ 1*Ability

           ## Label the error variances of the dependent variables
           Job_knowledge ~~ VarE_J*Job_knowledge
           Work_sample ~~ VarE_W*Work_sample
           Supervisor ~~ VarE_S*Supervisor"

RAM2 <- lavaan2RAM(model2, obs.variables=c("Ability","Job_knowledge",
                   "Work_sample","Supervisor"))

## The model is okay.     
checkRAM(Amatrix=RAM2$A, Smatrix=RAM2$S)   
}
}
\keyword{utilities}

