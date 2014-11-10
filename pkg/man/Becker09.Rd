\name{Becker09}
\alias{Becker09}
\docType{data}
\title{Ten Studies of Correlation Matrices used by Becker (2009) 
}
\description{
This data set includes ten studies of correlation matrices used by
Becker (2009).
}
\usage{data("Becker09")}

\details{
A list of data with the following structure:
\describe{
  \item{data}{A list of 4x4 correlation matrices. The
	variables are \emph{Performance}, \emph{Cognitive}, \emph{Somatic},
	and \emph{Self confidence}}
  \item{n}{A vector of sample sizes}
  \item{Type_of_sport}{Samples based on \emph{Individual} or \emph{Team}}
}
}

\source{
Craft, L. L., Magyar, T. M., Becker, B. J., & Feltz, D. L. (2003). The relationship between the Competitive State Anxiety Inventory-2 and sport performance: a meta-analysis. \emph{Journal of Sport and Exercise Psychology}, \bold{25(1)}, 44-65.
}
\references{
Becker, B. J. (2009). Model-based meta-analysis. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta-analysis} (2nd ed., pp. 377-395). New York: Russell Sage Foundation.
}
\examples{
data(Becker09)

#### Fixed-effects model
## First stage analysis
fixed1 <- tssem1(Becker09$data, Becker09$n, method="FEM")
summary(fixed1)

## Prepare a regression model using create.mxMatrix()
A1 <- create.mxMatrix(c(0, "0.1*Cog2Per", "0.1*SO2Per", "0.1*SC2Per",
                        0, 0, 0, 0,
                        0, 0, 0, 0,
                        0, "0.1*Cog2SC", "0.1*SO2SC",0),
                      type="Full", byrow=TRUE, ncol=4, nrow=4)

S1 <- create.mxMatrix(c("0.1*var_Per",
                        0, 1,
                        0, "0.1*cor", 1,
                        0, 0, 0, "0.1*var_SC"), byrow=TRUE, type="Symm")

## Second stage analysis
fixed2 <- tssem2(fixed1, Amatrix=A1, Smatrix=S1, diag.constraints=TRUE,
                 intervals.type="LB", model.name="TSSEM2 Becker09",
                 mx.algebras=list( Cog=mxAlgebra(Cog2SC*SC2Per, name="Cog"), 
                                   SO=mxAlgebra(SO2SC*SC2Per, name="SO"),
                                   Cog_SO=mxAlgebra(Cog2SC*SC2Per+SO2SC*SC2Per,
                                   name="Cog_SO")) )

## Rerun to clear the condition code
fixed2 <- rerun(fixed2)

summary(fixed2)


#### Fixed-effects model: with type of sport as cluster
## First stage analysis
cluster1 <- tssem1(Becker09$data, Becker09$n, method="FEM",
                   cluster=Becker09$Type_of_sport)
summary(cluster1)

## Second stage analysis  
cluster2 <- tssem2(cluster1, Amatrix=A1, Smatrix=S1, diag.constraints=TRUE,
                 intervals.type="LB", model.name="TSSEM2 Becker09",
                 mx.algebras=list( Cog=mxAlgebra(Cog2SC*SC2Per, name="Cog"), 
                                   SO=mxAlgebra(SO2SC*SC2Per, name="SO"),
                                   Cog_SO=mxAlgebra(Cog2SC*SC2Per+SO2SC*SC2Per,
                                   name="Cog_SO")) )
summary(cluster2)

#### Random-effects model
## First stage analysis
random1 <- tssem1(Becker09$data, Becker09$n, method="REM", RE.type="Diag")
summary(random1)

## Second stage analysis
random2 <- tssem2(random1, Amatrix=A1, Smatrix=S1, diag.constraints=TRUE,
                  intervals.type="LB", model.name="TSSEM2 Becker09",
                  mx.algebras=list( Cog=mxAlgebra(Cog2SC*SC2Per, name="Cog"), 
                                    SO=mxAlgebra(SO2SC*SC2Per, name="SO"),
                                    Cog_SO=mxAlgebra(Cog2SC*SC2Per+SO2SC*SC2Per,
                                    name="Cog_SO")) )
summary(random2)
}
\keyword{datasets}
