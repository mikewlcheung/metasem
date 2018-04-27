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

## Multiple-endpoint studies
Gleser94$MES
}
}
\keyword{datasets}
