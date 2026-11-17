# FIXME: ensure dimnames in Cov iputs?


#' Convert \code{metaSEM} objects into \code{semPlotModel} objects for plotting
#' 
#' It converts objects in class \code{wls} into objects of class
#' \code{semPlotModel}.
#' 
#' It uses the \code{ramModel()} to do the conversion.
#' 
#' @param object An object of class \code{wls} returned from \code{wls()} or
#' \code{tssem2()}.
#' @param manNames A character vector of the manifest names. The program will
#' try to get it from the \code{object} if it is not given.
#' @param latNames A character vector of the latent names. The program will
#' create it by using "L1", "L2", etc if it is not given.
#' @param labels Either \code{labels} (default if missing) or \code{RAM}. If
#' \code{labels}, the labels of the parameters are used in plotting. If
#' \code{RAM}, the RAM notations are used in plotting.
#' @param \dots Further arguments to be passed to
#' \code{\link[semPlot]{ramModel}}
#' @return A "semPlotModel" object.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[semPlot]{ramModel}}, \code{\link[metaSEM]{Becker92}},
#' \code{\link[metaSEM]{Becker09}}, \code{\link[metaSEM]{Digman97}},
#' \code{\link[metaSEM]{Hunter83}}
#' @keywords methods tssem
meta2semPlot <- function(object, manNames=NULL, latNames=NULL, labels=c("labels", "RAM"), ...) {

    if (!requireNamespace("semPlot", quietly=TRUE))    
        stop("\"semPlot\" package is required for this function.")
    
    if (inherits(object, "wls")) {
    ## No visible binding for global variable Note in R CMD check    
    ## A <- mxEval(Amatrix, object$mx.fit)
    ## S <- mxEval(Smatrix, object$mx.fit)
    ## F <- mxEval(Fmatrix, object$mx.fit)
       A <- unname(object$mx.fit@matrices$Amatrix$values)
       S <- object$mx.fit@matrices$Smatrix$values
       ## S can be calculated when diag.constraints=FALSE
       if (is.null(S)) {
           S <- object$mx.fit$algebras$Smatrix$result
       }
       S <- unname(S)

       F <- object$mx.fit@matrices$Fmatrix$values
       ## Try to get the names of the latent variables
       selector <- all(is.finite(F)) && all(F %in% c(0, 1)) &&
                   all(rowSums(F) == 1) && all(colSums(F) <= 1)
       observed <- if (selector) colSums(F) == 1 else
                   dimnames(F)[[2]] %in% dimnames(F)[[1]]
       lat_names <- colnames(F)[!observed]

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
        
       if (is.null(latNames)) {
           if (!is.null(lat_names)) {
               latNames <- lat_names
           } else {
               no.latent <- ncol(F) - nrow(F)
               if (no.latent > 0) latNames <- paste("L", seq(1,no.latent), sep="")
           }    
       }

       if (selector && !("Names" %in% names(list(...)))) {
           node.names <- character(ncol(F))
           node.names[observed] <- manNames[max.col(t(F[, observed, drop=FALSE]))]
           node.names[!observed] <- latNames
           out <- semPlot::ramModel(A=A, S=S, F=unname(F), manNames=manNames,
                                    latNames=latNames, Names=node.names,
                                    ObsCovs=ObsCovs, ImpCovs=ImpCovs, ...)
       } else {
           out <- semPlot::ramModel(A=A, S=S, F=unname(F), manNames=manNames,
                                    latNames=latNames, ObsCovs=ObsCovs,
                                    ImpCovs=ImpCovs, ...)
       }
    
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

  if (!is.element("wls", class(x)))
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

  if (!inherits(x, c("osmasem", "osmasem2", "osmasem3L")))
      stop("\"x\" must be an object of either class \"osmasem\", \"osmasem2\", or \"osmasem3L\".")

  at.means <- .osmasem.plot.at.means(x$mx.fit)
  get.matrix <- function(name, source=NULL)
    .osmasem.plot.matrix(name, at.means, source)
  A <- get.matrix("Amatrix", x$RAM$A)
  S <- get.matrix("Smatrix", x$RAM$S)

  F <- x$mx.fit$Fmatrix$values
  Id <- diag(nrow(S))

  ## Approximate sample covariance matrix
  ## The means in osmasem are correlation matrix
  if (inherits(x, "osmasem2")) {
    ObsCovs <- vec2symMat(apply(x$data$data[, x$data$ylabels, drop=FALSE], 2, mean, na.rm=TRUE),
                          diag=!x$cor.analysis)       
    ImpCovs <- get.matrix("expSigma")
    if (is.null(manNames)) manNames <- x$data$obslabels
  } else {
    ## A labels slot is required because osmasem allows subset.
    ## labels can be shorter than those in data$labels
    ObsCovs <- vec2symMat(apply(x$data$data[, x$labels$ylabels, drop=FALSE], 2, mean,
                                na.rm=TRUE), diag=FALSE)        
    ImpCovs <- get.matrix("impliedR")
    if (is.null(manNames)) manNames <- x$labels$obslabels
  }
  ## If manNames is still NULL
  if (is.null(manNames)) manNames <- paste("X", seq(1, nrow(F), 1), sep="")

  dimnames(ImpCovs) <- dimnames(ObsCovs) <- list(manNames, manNames)

  ## Map each observed RAM column to the corresponding row of F. This also
  ## handles data names and custom names that differ from the RAM names.
  selector <- all(is.finite(F)) && all(F %in% c(0, 1)) &&
              all(rowSums(F) == 1) && all(colSums(F) <= 1)
  observed <- if (selector) colSums(F) == 1 else
              colnames(F) %in% rownames(F)
  lat_names <- colnames(F)[!observed]
  if (is.null(latNames)) {
    if (!is.null(lat_names)) latNames <- lat_names
    else if (ncol(F) > nrow(F)) latNames <- paste0("L", seq_len(ncol(F)-nrow(F)))
  }
  args <- list(A=A, S=S, F=F, manNames=manNames, latNames=latNames,
               ObsCovs=ObsCovs, ImpCovs=ImpCovs)
  if (selector) {
    if (length(manNames) != sum(observed) || length(latNames) != sum(!observed))
      stop("Lengths of 'manNames' and 'latNames' must match the observed and latent variables.")
    node.names <- character(ncol(F))
    node.names[observed] <- manNames[max.col(t(F[, observed, drop=FALSE]))]
    node.names[!observed] <- latNames
    args$Names <- node.names
  }
  if (inherits(x, "osmasem2") && x$mean.analysis) {
    M <- get.matrix("Mmatrix", x$RAM$M)
    colnames(M) <- colnames(F)
    args$M <- M
  }
  out <- do.call(semPlot::ramModel, args)

  if (inherits(x, c("osmasem", "osmasem3L"))) {
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


## Use the fitted rows, including any subset, as the reference population.
.osmasem.plot.at.means <- function(model) {
  if (is.null(model@data) || model@data$type != "raw") return(model)
  variables <- unique(unlist(lapply(model@matrices,function(m) c(m$labels))))
  variables <- sub("^data\\.","",grep("^data\\.",variables,value=TRUE))
  observed <- model@data@observed
  variables <- intersect(variables,names(observed))
  if (!length(variables)) return(model)
  reference <- observed[1,,drop=FALSE]
  for (v in variables) reference[[v]] <- mean(observed[[v]],na.rm=TRUE)
  model@data <- mxData(reference,type="raw",primaryKey=model@data@primaryKey)
  for (nm in names(model@matrices)) {
    m <- model@matrices[[nm]]
    for (v in variables) {
      index <- !is.na(m$labels) & m$labels == paste0("data.",v)
      m$values[index] <- reference[[v]]
    }
    model@matrices[[nm]] <- m
  }
  model
}

.osmasem.plot.matrix <- function(name, model, source=NULL) {
  object <- model[[name]]
  if (is.null(object)) return(NULL)
  value <- if (inherits(object,"MxMatrix")) object$values else
    mxEvalByName(name,model,compute=TRUE)
  if (is.null(dimnames(value))) {
    fallback <- dimnames(source)
    if (is.null(fallback) && !is.null(model$Fmatrix)) {
      all.names <- colnames(model$Fmatrix$values)
      observed <- rownames(model$Fmatrix$values)
      fallback <- switch(name,
        Amatrix=list(all.names,all.names), Smatrix=list(all.names,all.names),
        Mmatrix=list(NULL,all.names), expSigma=list(observed,observed),
        impliedR=list(observed,observed), NULL)
    }
    dimnames(value) <- fallback
  }
  value
}
