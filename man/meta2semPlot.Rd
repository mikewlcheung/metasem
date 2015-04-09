\name{meta2semPlot}
\alias{meta2semPlot}
\title{
Convert \code{metaSEM} objects into \code{semPlotModel} objects for plotting
}
\description{
It converts objects in class \code{wls} into objects of class \code{semPlotModel}.
}
\usage{
meta2semPlot(object, manNames = NULL, latNames = NULL, labels = c("labels", "RAM"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{wls} returned from \code{wls()}
	or \code{tssem2()}.
}
  \item{manNames}{A character vector of the manifest names. The program
	will try to get it from the \code{object} if it is not given. 
}
  \item{latNames}{A character vector of the latent names. The program
	will create it by using "L1", "L2", etc if it is not given. 
}
  \item{labels}{Either \code{labels} (default if missing) or
	\code{RAM}. If \code{labels}, the labels of the parameters are used
	in plotting. If \code{RAM}, the RAM notations are used in plotting.
}
  \item{\dots}{Further arguments to be passed to \code{\link[semPlot]{ramModel}}
}
}
\details{
It uses the \code{ramModel()} to do the conversion.
}
\value{
  A "semPlotModel" object.
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}

\seealso{\code{\link[semPlot]{ramModel}},
  \code{\link[metaSEM]{Becker92}}, \code{\link[metaSEM]{Becker09}},
  \code{\link[metaSEM]{Digman97}}, \code{\link[metaSEM]{Hunter83}}
}

\keyword{ methods }
\keyword{ tssem }% __ONLY ONE__ keyword per line
