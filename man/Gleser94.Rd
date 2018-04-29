\name{Gleser94}
\alias{Gleser94}
\docType{data}
\title{Two Datasets from Gleser and Olkin (1994)
}
\description{It includes two datasets in multiple-treatment studies and
  multiple-endpoint studies reported by Gleser and Olkin (1994).
}
\usage{data("Gleser94")}
\format{
  A list of two data frames.
  \describe{
    \item{\code{MTS}}{A data frame of multiple-treatment studies.}
    \item{\code{MES}}{A data frame of multiple-endpoint studies.}	
  }
}
\source{Gleser, L. J., & Olkin, I. (1994). Stochastically dependent effect sizes. In H. Cooper & L. V. Hedges (Eds.), The handbook of research synthesis. (pp. 339-355). New York: Russell Sage Foundation.
}

\seealso{\code{\link[metaSEM]{smdMTS}}, \code{\link[metaSEM]{smdMES}}
}

\examples{
\dontrun{
data(Gleser94)

## Multiple-treatment studies
Gleser94$MTS

means <- paste0("Mean.", c("C", "E1", "E2", "E3", "E4", "E5"))
sds <-  paste0("SD.", c("C", "E1", "E2", "E3", "E4", "E5"))
ns <- paste0("N.", c("C", "E1", "E2", "E3", "E4", "E5"))

## Assuming homogeneity of variances
t(apply(Gleser94$MTS, MARGIN=1,
        function(x)
            smdMTS(m=x[c("Mean.C", "Mean.E1", "Mean.E2", "Mean.E3", "Mean.E4", "Mean.E5")],
                   v=x[c("SD.C", "SD.E1", "SD.E2", "SD.E3", "SD.E4", "SD.E5")]^2,
                   n=x[c("N.C", "N.E1", "N.E2", "N.E3", "N.E4", "N.E5")],
                   homogeneity="variance", list.output=FALSE)))
                          
## Multiple-endpoint studies
Gleser94$MES

## Calculate the sampling variances and covariance
Gleser94$MES$Uncoached.V11 <- with(Gleser94$MES, SD.Uncoached.Math^2)
Gleser94$MES$Uncoached.V21 <- with(Gleser94$MES, SD.Uncoached.Math*Cor.Math.Verbal*SD.Uncoached.Verbal)
Gleser94$MES$Uncoached.V22 <- with(Gleser94$MES, SD.Uncoached.Verbal^2)

Gleser94$MES$Coached.V11 <- with(Gleser94$MES, SD.Coached.Math^2)
Gleser94$MES$Coached.V21 <- with(Gleser94$MES, SD.Coached.Math*Cor.Math.Verbal*SD.Coached.Verbal)
Gleser94$MES$Coached.V22 <- with(Gleser94$MES, SD.Coached.Verbal^2)


## Assuming homogeneity of covariance matrices
t(apply(Gleser94$MES, MARGIN=1,
        function(x)
            smdMES(m1=x[c("Mean.Uncoached.Math", "Mean.Uncoached.Verbal")],
                   m2=x[c("Mean.Coached.Math", "Mean.Coached.Verbal")],
                   V1=vec2symMat(x[c("Uncoached.V11", "Uncoached.V21", "Uncoached.V22")]),
                   V2=vec2symMat(x[c("Coached.V11", "Coached.V21", "Coached.V22")]),
                   n1=x["N.Uncoached"],
                   n2=x["N.Coached"],
                   homogeneity="covariance", list.output=FALSE)))
}
}
\keyword{datasets}
