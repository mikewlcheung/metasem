# FIXME: ensure dimnames in Cov iputs?
meta2semPlot <- function(object, manNames=NULL, latNames=NULL, labels=c("labels", "RAM"), ...) {

    if (!requireNamespace("semPlot", quietly=TRUE))    
        stop("\"semPlot\" package is required for this function.")
    
    if (inherits(object, "wls")) {
    ## No visible binding for global variable Note in R CMD check    
    ## A <- mxEval(Amatrix, object$mx.fit)
    ## S <- mxEval(Smatrix, object$mx.fit)
    ## F <- mxEval(Fmatrix, object$mx.fit)
       A <- object$mx.fit@matrices$Amatrix$values
       S <- object$mx.fit@matrices$Smatrix$values
       ## S can be calculated when diag.constraints=FALSE
       if (is.null(S)) {
           S <- object$mx.fit$algebras$Smatrix$result
           ## Add the dimensions names. Not necessary, but nicer
           ## dimnames(S) <- dimnames(A)
       }
       F <- object$mx.fit@matrices$Fmatrix$values
       Id <- diag(nrow(S))
       ObsCovs <- object$Cov
    ## ImpCovs <- mxEval(Fmatrix%*%solve(Id-Amatrix)%*%S%*%solve(t(Id-Amatrix))%*%t(Fmatrix), 
    ##                   object$mx.fit)
    ## ImpCovs <- F%*%solve(Id-A)%*%S%*%solve(t(Id-A))%*%t(F)
        
       ImpCovs <- object$mx.fit@algebras$impliedS$result    
        
       if (is.null(manNames)) {
               manNames <- colnames(ObsCovs)
       }    
       ## If colnames(ObsCovs) is still NULL
       if (is.null(manNames)) {
           manNames <- paste("X", seq(1, nrow(F), 1), sep="")
       }

       dimnames(ImpCovs) <- dimnames(ObsCovs) <- list(manNames, manNames)

       ## Try to get the names of the latent variables
       lat_index <- dimnames(F)[[2]] %in% dimnames(F)[[1]]
       lat_names <- dimnames(F)[[2]][!lat_index]
        
       if (is.null(latNames)) {
           if (!is.null(lat_names)) {
               latNames <- lat_names
           } else {
               no.latent <- ncol(F) - nrow(F)
               if (no.latent > 0) latNames <- paste("L", seq(1,no.latent), sep="")
           }    
       }
    
       out <- semPlot::ramModel(A=A, S=S, F=F, manNames=manNames, latNames=latNames, ObsCovs=ObsCovs, 
                                ImpCovs=ImpCovs, ...)
    
       labels <- match.arg(labels)    
       if (labels=="labels") {
         A.labels <- object$mx.fit@matrices$Amatrix$labels
         S.labels <- object$mx.fit@matrices$Smatrix$labels
         labels <- c(c(A.labels), c(S.labels))
         row.pars <- as.numeric(row.names(out@Pars))
         out@Pars$label <- labels[row.pars]
         ## replace NA labels with space
         na.index <- is.na(out@Pars$label)
         ##out@Pars$label[na.index] <- sprintf("%.2f", out@Pars$est[na.index])
         out@Pars$label[na.index] <- ""
       }   
  } else {
    stop("\"object\" must be an object of class \"wls\".")
    
  }
    
  out
}

plot.character <- function(x, fixed.x=FALSE, nCharNodes=0, nCharEdges=0,
                           layout=c("tree", "circle", "spring", "tree2", "circle2"),
                           sizeMan=8, sizeLat=8, edge.label.cex=1.3,
                           color="white", ...) {
  if (!requireNamespace("semPlot", quietly=TRUE))    
    stop("\"semPlot\" package is required for this function.")
  
  if (!is.element("character", class(x)))
    stop("\"x\" must be an object of class \"character\" of the \"lavaan\" model.")

  invisible( semPlot::semPaths(semPlot::semPlotModel(x, fixed.x=fixed.x),
                               nCharNodes=nCharNodes, nCharEdges=nCharEdges,
                               layout=match.arg(layout), sizeMan=sizeMan,
                               sizeLat=sizeLat, edge.label.cex=edge.label.cex,
                               color=color, ...) )
}

## How to handle multiple plots with cluster?
plot.wls <- function(x, manNames=NULL, latNames=NULL, labels=c("labels", "RAM"),
                     what="est", nCharNodes=0, nCharEdges=0,
                     layout=c("tree", "circle", "spring", "tree2", "circle2"),
                     sizeMan=8, sizeLat=8, edge.label.cex=1.3, color="white",
                     weighted=FALSE, ...) {
  if (!requireNamespace("semPlot", quietly=TRUE))    
    stop("\"semPlot\" package is required for this function.")

  if (!is.element(class(x), "wls"))
    stop("\"x\" must be an object of class \"wls\".")

  sem.plot <- meta2semPlot(object=x, manNames=manNames, latNames=latNames, 
                           labels=labels)

  invisible( semPlot::semPaths(sem.plot, what=what, nCharNodes=nCharNodes,
                               nCharEdges=nCharEdges, layout=match.arg(layout),
                               sizeMan=sizeMan, sizeLat=sizeLat,
                               edge.label.cex=edge.label.cex, color=color,
                               weighted=weighted, ...) )
}

plot.osmasem <- function(x, manNames=NULL, latNames=NULL, labels=c("labels", "RAM"),
                         what="est", nCharNodes=0, nCharEdges=0,
                         layout=c("tree", "circle", "spring", "tree2", "circle2"),
                         sizeMan=8, sizeLat=8, edge.label.cex=1.3, color="white",
                         weighted=FALSE, ...) {
  if (!requireNamespace("semPlot", quietly=TRUE))    
    stop("\"semPlot\" package is required for this function.")

  if (!is.element(class(x), c("osmasem", "osmasem2", "osmasem3L")))
      stop("\"x\" must be an object of either class \"osmasem\", \"osmasem2\", or \"osmasem3L\".")

  if (is.element(class(x), "osmasem2")) {
    A <- x$mx.fit$Amatrix$values
    S <- x$mx.fit$Smatrix$values
  } else {
    ## osmasem and osmasem3L: different from osmasem2 coz they replace error
    ## variances with other parameters 
    A <- x$mx.fit$Amatrix$result
    S <- x$mx.fit$Smatrix$result
    ## S can be calculated when diag.constraints=FALSE
    if (is.null(S)) {
      S <- x$mx.fit$algebras$Smatrix$result
    }
  }

  F <- x$mx.fit$Fmatrix$values
  Id <- diag(nrow(S))

  ## Approximate sample covariance matrix
  ## The means in osmasem are correlation matrix
  if (is.element(class(x), "osmasem2")) {    
    ObsCovs <- vec2symMat(apply(x$data$data[, x$data$ylabels], 2, mean, na.rm=TRUE),
                          diag=!x$cor.analysis)       
    ImpCovs <- x$mx.fit$expSigma$result
    if (is.null(manNames)) manNames <- x$data$obslabels
  } else {
    ## A labels slot is required because osmasem allows subset.
    ## labels can be shorter than those in data$labels
    ObsCovs <- vec2symMat(apply(x$data$data[, x$labels$ylabels], 2, mean,
                                na.rm=TRUE), diag=FALSE)        
    ImpCovs <- x$mx.fit$impliedR$result
    if (is.null(manNames)) manNames <- x$labels$obslabels
  }
  ## If manNames is still NULL
  if (is.null(manNames)) manNames <- paste("X", seq(1, nrow(F), 1), sep="")

  dimnames(ImpCovs) <- dimnames(ObsCovs) <- list(manNames, manNames)

  ## Try to get the names of the latent variables
  lat_index <- dimnames(F)[[2]] %in% dimnames(F)[[1]]
  lat_names <- dimnames(F)[[2]][!lat_index]

  if (is.null(latNames)) {
      if (!is.null(lat_names)) {
          latNames <- lat_names
      } else {
          no.latent <- ncol(F) - nrow(F)
          if (no.latent > 0) latNames <- paste("L", seq(1,no.latent), sep="")
      }    
  }

  if (is.element(class(x), "osmasem2")) {
    if (x$mean.analysis) {
      ## osmasem2 with mean structure
      M <- x$mx.fit$Mmatrix$values
      colnames(M) <- colnames(F)
      out <- semPlot::ramModel(A=A, S=S, F=F, M=M, manNames=manNames,
                               latNames=latNames, ObsCovs=ObsCovs, 
                               ImpCovs=ImpCovs)
    } else {
      ## osmasem2 without mean structure
      out <- semPlot::ramModel(A=A, S=S, F=F, manNames=manNames, latNames=latNames,
                               ObsCovs=ObsCovs, ImpCovs=ImpCovs)
    }
  } else {
    ## osmasem and osmasem3L
    out <- semPlot::ramModel(A=A, S=S, F=F, manNames=manNames, latNames=latNames,
                             ObsCovs=ObsCovs, ImpCovs=ImpCovs)
  }

  if (is.element(class(x), c("osmasem", "osmasem3L"))) {
    labels <- match.arg(labels)
    if (labels=="labels") {
      A.labels <- x$mx.fit$A0$labels
      S.labels <- x$mx.fit$S0$labels
      labels <- c(c(A.labels), c(S.labels))
      row.pars <- as.numeric(row.names(out@Pars))
      out@Pars$label <- labels[row.pars]
      ## replace NA labels with space
      na.index <- is.na(out@Pars$label)
      ##out@Pars$label[na.index] <- sprintf("%.2f", out@Pars$est[na.index])
      out@Pars$label[na.index] <- ""
    }
  }

  invisible( semPlot::semPaths(out, what=what, nCharNodes=nCharNodes,
                               nCharEdges=nCharEdges, layout=match.arg(layout),
                               sizeMan=sizeMan, sizeLat=sizeLat,
                               edge.label.cex=edge.label.cex, color=color,
                               weighted=weighted, ...) )
}

plot.osmasem2 <- function(x, manNames=NULL, latNames=NULL, labels=c("labels", "RAM"),
                          what="est", nCharNodes=0, nCharEdges=0,
                          layout=c("tree", "circle", "spring", "tree2", "circle2"),
                          sizeMan=8, sizeLat=8, edge.label.cex=1.3, color="white",
                          weighted=FALSE, ...) {
  plot.osmasem(x, manNames=manNames, latNames=latNames, labels=labels,
               what=what, nCharNodes=nCharNodes, nCharEdges=nCharEdges,
               layout=layout, sizeMan=sizeMan, sizeLat=sizeLat,
               edge.label.cex=edge.label.cex, color=color, weighted=weighted, ...)
}

## plot.osmasem3L <- function(x, manNames=NULL, latNames=NULL, labels=c("labels", "RAM"),
##                            what="est", nCharNodes=0, nCharEdges=0,
##                            layout=c("tree", "circle", "spring", "tree2", "circle2"),
##                            sizeMan=8, sizeLat=8, edge.label.cex=1.3, color="white",
##                            weighted=FALSE, ...) {
##     plot.osmasem(x, manNames=manNames, latNames=latNames, labels=labels,
##                  what=what, nCharNodes=nCharNodes, nCharEdges=nCharEdges,
##                  layout=layout, sizeMan=sizeMan, sizeLat=sizeLat,
##                  edge.label.cex=edge.label.cex, color=color, weighted=weighted, ...)
## }
