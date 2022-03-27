rCor3L <- function(Sigma, V.B, V.W, n, cluster, corr=TRUE, raw.data=FALSE,
                   nonPD.pop=c("replace", "nearPD", "accept"),
                   nonPD.sam=c("stop", "nearPD")) {
  
    nonPD.pop <- match.arg(nonPD.pop)
    nonPD.sam <- match.arg(nonPD.sam)

    if (sum(cluster) != length(n)) {
        stop("The length of 'n' must equal the sum of 'cluster'.\n")
    }
  
    ## Generate between-level population matrices   
    P.B <- rCorPop(Sigma=Sigma, V=V.B, k=length(cluster), corr=corr, nonPD.pop=nonPD.pop)

    ## Generate within-level population matrices
    P.W <- mapply(rCorPop, Sigma=P.B, k=cluster, 
                  MoreArgs=list(V=V.W, corr=corr, nonPD.pop=nonPD.pop),
                  SIMPLIFY=FALSE)
  
    names(P.B) <- names(P.W) <- paste0("Cluster", seq_along(cluster))  
  
    ## Add "_" in the names, which make it easier to read later
    P.W.tmp <- P.W
    names(P.W.tmp) <- paste0("Cluster", seq_along(cluster), "_") 
  
    ## Generate sample matrices
    R <- rCorSam(Sigma=unlist(P.W.tmp, recursive=FALSE), n=n, corr=corr, 
                 raw.data=raw.data, nonPD.sam=nonPD.sam)
  
    ## Labels for the clusters
    Cluster <- paste0("Cluster", rep(seq_along(cluster), times=cluster))
  
    out <- list(P.B=P.B, P.W=P.W, R=R, cluster=Cluster, n=n)
    attr(out, "Sigma") <- Sigma
    attr(out, "V.B") <- V.B
    attr(out, "V.W") <- V.W
    attr(out, "cluster") <- cluster
    class(out) <- "Cor3L"
    out    
}

summary.Cor3L <- function(object, ...) {
  if (!is.element("Cor3L", class(object)))
    stop("\"object\" must be an object of class \"Cor3L\".")  

  cluster <- attr(object, "cluster")
  
  sum.b <- summary(object$P.B)
  sum.w <- lapply(object$P.W, summary)
  
  ## Numbers of within studies. Should be the same as cluster.
  k.w <- sapply(sum.w, function(x) x$k)

  ## Empirical V.W
  V.W.emp <- lapply(sum.w, function(x) x$V_Samp)
  V.W.emp <- Reduce("+", Map("*", k.w, V.W.emp))/sum(k.w)
  nonPD.pop.W <- sum.w[[1]]$nonPD.pop
  nonPD.count.W <- sum(sapply(sum.w, function(x) x$nonPD.count))
 
  out <- list(Sigma=attr(object, "Sigma"),
              V.B = attr(object, "V.B"),
              V.W = attr(object, "V.W"),
              cluster = cluster,
              Sigma.emp = sum.b$R,
              V.B.emp = sum.b$V_Samp,
              nonPD.pop.B = sum.b$nonPD.pop,
              nonPD.count.B = sum.b$nonPD.count,
              V.W.emp = V.W.emp,
              nonPD.pop.W = nonPD.pop.W,
              nonPD.count.W = nonPD.count.W)
  class(out) <- "summary.Cor3L"
  out
}

print.summary.Cor3L <- function(x, ...) {
  if (!is.element("summary.Cor3L", class(x)))
    stop("\"x\" must be an object of class \"summary.Cor3L\".")
  
  cat("Population Sigma:\n")
  print(x$Sigma)
  cat("\nCluster sizes:\n")
  print(x$cluster)
  
  cat("\nPopulation V (between):\n")
  print(x$V.B)
  cat("\nEmpirical V (between):\n")
  print(x$V.B.emp)
  cat("\nMethod to handle non-positive definite matrices (between):", x$nonPD.pop.B)
  cat("\nNumber of samples (between):", length(x$cluster))
  cat("\nCount of non-positive definite matrices (between):", x$nonPD.count.B, "\n")

  cat("\nPopulation V (within):\n")
  print(x$V.W)
  cat("\nEmpirical V (within):\n")
  print(x$V.W.emp)
  cat("\nMethod to handle non-positive definite matrices (within):", x$nonPD.pop.W)
  cat("\nNumber of samples (within):", sum(x$cluster))
  cat("\nCount of non-positive definite matrices (within):", x$nonPD.count.W, "\n")
}  
              
