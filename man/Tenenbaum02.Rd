\name{Tenenbaum02}
\alias{Tenenbaum02}
\docType{data}
\title{Correlation coefficients reported by Tenenbaum and Leaper (2002)}
\description{
  Fourty-eight studies reported by Tenenbaum and Leaper (2002, Table 1).
}
\usage{data(Tenenbaum02)}
\details{
    The variables are:
  \describe{
  \item{Authors}{Authors of the study}
  \item{Year}{Year of publication}	
  \item{N}{Sample size}
  \item{r}{Correlation between parents' gender schemas and their
	offspring's gender-related cognitions.}
  \item{v}{Sampling variance of r}
  \item{Publication_source}{Publication source: 1="top-tier journal",
	2="second-tier journal or book chapter", 3="dissertation", 4="other unpublished study"}
  \item{Author_gender}{Gender of the first author: "W"="woman", "M"="man"}
  \item{Parent_type}{Parent type: "M"="mother", "F"="father", "MF"="mother and father"}
  \item{Parent_predictor}{Parent predictor: "S"="self gender schema", "A"="gender attitudes about others"}
  \item{Offspring_age}{Offspring age (months)}
  \item{Offspring_type}{Offspring type: "D"="daughter",
	"S"="son", "DS"="daughter and son"}
  \item{Offspring_outcome}{Offspring outcome:
	"S"="gender schema for self", "A"="gender attitudes toward others",
	"I"="gender-related interests and preferences", "W"="work-related attitudes"}  
}

}
\source{
  Tenenbaum, H. R., & Leaper, C. (2002). Are parents' gender schemas related to their children's gender-related cognitions? A meta-analysis. \emph{Developmental Psychology}, \emph{38}(4), 615-630. https://doi.org/10.1037/0012-1649.38.4.615
}

\examples{
\dontrun{
data(Tenenbaum02)
}
}
\keyword{datasets}
