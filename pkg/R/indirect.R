standIndirect <- function(x, n) {
  if (is.list(x)) {
    t(mapply(standIndirect, x, n=n))
  } else {
  dimnames(x) <- list(c("y","m","x"),c("y","m","x"))
  ## starting values for standard deviations
  stvalues <- sqrt(diag(x))
  myMod <- mxModel("Standardized indirect effect", type="RAM", mxData(observed=x, type="cov", numObs=n),
                   manifestVars=c("y","m","x"),
                   ## Im: imaginary variables
                   latentVars=c("yLat","mLat","xLat","xIm1","mIm2","xIm2","varIm1","varIm2","covIm"),
                   mxPath(from="xLat", to="x", arrows=1, free=TRUE, values=stvalues[3], labels="xsd"),
                   mxPath(from="mLat", to="m", arrows=1, free=TRUE, values=stvalues[2], labels="msd"),                   
                   mxPath(from="yLat", to="y", arrows=1, free=TRUE, values=stvalues[1], labels="ysd"),
                   mxPath(from="xLat", arrows=2, free=FALSE, values=1),
                   mxPath(from="mLat", arrows=2, free=FALSE, values=0),
                   mxPath(from="yLat", arrows=2, free=FALSE, values=0),                   
                   mxPath(from="xIm1", arrows=2, free=FALSE, values=-1),
                   mxPath(from="xIm2", arrows=2, free=FALSE, values=-1),
                   mxPath(from="mIm2", arrows=2, free=FALSE, values=-1),
                   mxPath(from="varIm1", arrows=2, free=FALSE, values=1),
                   mxPath(from="varIm2", arrows=2, free=FALSE, values=1),
                   mxPath(from="covIm", arrows=2, free=FALSE, values=0),
                   mxPath(from="x", arrows=2, free=FALSE, values=0),
                   mxPath(from="m", arrows=2, free=FALSE, values=0),
                   mxPath(from="y", arrows=2, free=FALSE, values=0),
                   mxPath(from="xLat", to="mLat", arrows=1, free=TRUE, values=0.2, labels="a"),
                   mxPath(from="mLat", to="yLat", arrows=1, free=TRUE, values=0.2, labels="b"),
                   mxPath(from="xLat", to="yLat", arrows=1, free=TRUE, values=0.2, labels="c"),
                   mxPath(from=c("varIm1","xIm1"), to="mLat", free=c(FALSE, TRUE), values=c(1,0.2), labels=c(NA,"a")),
                   mxPath(from=c("varIm2","xIm2","mIm2"), to="yLat", free=c(FALSE, TRUE, TRUE),
                          values=c(1,0.2,0.2), labels=c(NA,"c","b")),
                   mxPath(from="xIm2", to="covIm", arrows=2, free=TRUE, values=0.2, labels="a"),
                   mxPath(from="covIm", to="mIm2", arrows=2, free=FALSE, values=-1))
  
  my.fit <- mxRun(myMod, silent=TRUE)
  my.summary <- summary(my.fit)
  a <- mxEval(a, my.fit)
  b <- mxEval(b, my.fit)
  a.se <- my.summary$parameters[6,6]
  b.se <- my.summary$parameters[3,6]
  indirect <- a*b
  indirect.var <- a^2*b.se^2+b^2*a.se^2
  acovS <- tryCatch( 2*solve(my.fit@output$calculatedHessian[c("a","b","c"), c("a","b","c")]), error = function(e) e)
  if (inherits(acovS, "error")) {
    cat("Asymptotic covariance matrix of the estimates is not positive definite.\n")
    stop(print(acovS))
  }  

  c(ind.eff=indirect, dir.eff=mxEval(c, my.fit), ind.var=indirect.var, ind.dir.cov=a*acovS[3,2],
    dir.var=my.summary$parameters[5,6]^2)
  }
}
