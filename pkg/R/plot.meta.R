library(metaSEM)
library(ellipse)
library(MASS)
library(bayesm)
library(metafor)

## effect.sizes=x,y
plot.meta <- function(x, effect.sizes=seq(1, x$no.y), add.margin=0.1, interval=0.95,
                      main=paste("Effect Sizes and ",interval*100,"% Confidence Ellipses",sep=""),
                      axis.label=paste("Effect size ", effect.sizes, sep=""),
                      
                      study.col="black", study.pch=19, min.cex=0.8, weight.plot=FALSE,
                      study.ellipse.plot=TRUE, study.ellipse.col="black",
                      study.ellipse.lty=2, study.ellipse.lwd=0.5,

                      estimate.col="black", estimate.pch=18, estimate.cex=2,
                      estimate.ellipse.plot=TRUE, estimate.ellipse.col="black",
                      estimate.ellipse.lty=1, estimate.ellipse.lwd=2,

                      randeff.ellipse.plot=TRUE, randeff.ellipse.col="black",
                      randeff.ellipse.lty=1, randeff.ellipse.lwd=2,

                      univariate.plot=TRUE, univariate.lines.col="gray",
                      univariate.lines.lty=3, univariate.lines.lwd=1,
                      univariate.polygon.col="gray",
                      univariate.arrows.col="black", univariate.arrows.lwd=2, diag.panel=FALSE, ...) {

  if (!require(ellipse))
    stop("You need to install the \"ellipse\" package before you can use the plot.meta function.\n")
  if (!is.element("meta", class(x)))
    stop("\"object\" must be an x of class \"meta\".")
  
  no.y <- x$no.y
  if (x$no.y==1) stop("There must be at least TWO effect sizes.\n")
  if (x$no.x!=0) warning("There are predictors in the model.\nThe plot is based on the intercepts.\n")
  
  if (length(effect.sizes)>2) {
   
    if (diag.panel==FALSE) {
      mat <- vec2symMat(seq(1, no.y*(no.y-1)/2))
      mat[upper.tri(mat)] <- 0
      layout(mat, respect=TRUE)
      for (i in 1:(length(effect.sizes)-1))
        for (j in (i+1):length(effect.sizes)) {
          my.effects <- effect.sizes[c(i,j)]
          plot.meta(x, effect.sizes=my.effects,
                    add.margin=add.margin, interval=interval, main=main, axis.label=axis.label[c(i,j)],                    
                    study.col=study.col, study.pch=study.pch, min.cex=min.cex, weight.plot=weight.plot,
                    study.ellipse.plot=study.ellipse.plot, study.ellipse.col=study.ellipse.col,
                    study.ellipse.lty=study.ellipse.lty, study.ellipse.lwd=study.ellipse.lwd,
                    estimate.col=estimate.col, estimate.pch=estimate.pch, estimate.cex=estimate.cex,
                    estimate.ellipse.plot=estimate.ellipse.plot, estimate.ellipse.col=estimate.ellipse.col,
                    estimate.ellipse.lty=estimate.ellipse.lty, estimate.ellipse.lwd=estimate.ellipse.lwd,
                    randeff.ellipse.plot=randeff.ellipse.plot, randeff.ellipse.col=randeff.ellipse.col,
                    randeff.ellipse.lty=randeff.ellipse.lty, randeff.ellipse.lwd=randeff.ellipse.lwd,
                    univariate.plot=univariate.plot, univariate.lines.col=univariate.lines.col,
                    univariate.lines.lty=univariate.lines.lty, univariate.lines.lwd=univariate.lines.lwd,
                    univariate.polygon.col=univariate.polygon.col, univariate.arrows.col=univariate.arrows.col,
                    univariate.arrows.lwd=univariate.arrows.lwd, diag.panel=FALSE, ...)
        }              
    } else {
      mat <- vec2symMat(seq(1, no.y*(no.y-1)/2), diag=FALSE)
      mat[upper.tri(mat)] <- 0
      diag(mat) <- seq(no.y*(no.y-1)/2+1, no.y*(no.y+1)/2) 
      layout(mat, respect=TRUE)
      for (i in 1:(length(effect.sizes)-1))
        for (j in (i+1):length(effect.sizes)) {
          my.effects <- effect.sizes[c(i,j)]
          plot.meta(x, effect.sizes=my.effects,
                    add.margin=add.margin, interval=interval, main=main, axis.label=axis.label[c(i,j)],                    
                    study.col=study.col, study.pch=study.pch, min.cex=min.cex, weight.plot=weight.plot,
                    study.ellipse.plot=study.ellipse.plot, study.ellipse.col=study.ellipse.col,
                    study.ellipse.lty=study.ellipse.lty, study.ellipse.lwd=study.ellipse.lwd,
                    estimate.col=estimate.col, estimate.pch=estimate.pch, estimate.cex=estimate.cex,
                    estimate.ellipse.plot=estimate.ellipse.plot, estimate.ellipse.col=estimate.ellipse.col,
                    estimate.ellipse.lty=estimate.ellipse.lty, estimate.ellipse.lwd=estimate.ellipse.lwd,
                    randeff.ellipse.plot=randeff.ellipse.plot, randeff.ellipse.col=randeff.ellipse.col,
                    randeff.ellipse.lty=randeff.ellipse.lty, randeff.ellipse.lwd=randeff.ellipse.lwd,
                    univariate.plot=univariate.plot, univariate.lines.col=univariate.lines.col,
                    univariate.lines.lty=univariate.lines.lty, univariate.lines.lwd=univariate.lines.lwd,
                    univariate.polygon.col=univariate.polygon.col, univariate.arrows.col=univariate.arrows.col,
                    univariate.arrows.lwd=univariate.arrows.lwd, diag.panel=FALSE, ...)
        }       
    }
  # if (length(effect.sizes)>2)   
  } else {
    ## Only two effect sizes
    if ( (length(effect.sizes)==2)&(diag.panel==TRUE) ) {
      layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(1,1), c(1,1), respect=TRUE)
    }
      
    my.effects <- paste("Intercept", effect.sizes, sep="")
    ES <- coef(x)[my.effects]
    ACov <- vcov(x)[my.effects, my.effects]
    # Fixme: diff nos of rand.effects or fixed-effects only
    my.rand <- vech(outer(seq(1, length(my.effects)), seq(1, length(my.effects)),
                          function(x,y) { paste("Tau2_",x,"_",y,sep="")}))
    RE <- vec2symMat(coef(x)[my.rand])
    dimnames(RE) <- list(my.effects, my.effects)

  # y, x, var(y), cov(x,y), var(x)
    ## Index for positions of the conditional sampling variances
    rand.ind <- vec2symMat(seq(1, no.y*(no.y+1)/2))+ no.y
    ellipse.pt <- lapply(split(x$data, 1:nrow(x$data)),
                         function(x) {rand.val <- c(x[rand.ind[effect.sizes[1],effect.sizes[1]]],
                                                    x[rand.ind[effect.sizes[1],effect.sizes[2]]],
                                                    x[rand.ind[effect.sizes[2],effect.sizes[2]]])
                                      ellipse( vec2symMat(rand.val),
                                               centre=c(x[effect.sizes[1]],x[effect.sizes[2]]))})
    if (weight.plot==TRUE) {
      weigh <- sapply(split(x$data, 1:nrow(x$data)),                   
                      function(x) {rand.val <- c(x[rand.ind[effect.sizes[1],effect.sizes[1]]],
                                                 x[rand.ind[effect.sizes[1],effect.sizes[2]]],
                                                 x[rand.ind[effect.sizes[2],effect.sizes[2]]])
                                   sqrt(1/det( vec2symMat(rand.val) ))})
      study.cex <- weigh*min.cex/min(weigh)
    } else {
      study.cex <- min.cex
    }
    xlim <- sapply(ellipse.pt, function (x) { x <- data.frame(x); x$x })
    xlim <- c(min(xlim), max(xlim)) + c(-add.margin, 0)
    ylim <- sapply(ellipse.pt, function (x) { x <- data.frame(x); x$y })
    ylim <- c(min(ylim), max(ylim)) + c(-add.margin, 0)

    complete <- with(x, complete.cases(data[, effect.sizes[1]],data[, effect.sizes[2]]))
    my.x <- x$data[complete, effect.sizes[1]]
    my.y <- x$data[complete, effect.sizes[2]]
    plot(my.y~my.x, xlim=xlim, ylim=ylim, col=study.col, pch=study.pch, cex=study.cex,
         xlab=axis.label[1], ylab=axis.label[2], main=main, ...)
    if (study.ellipse.plot==TRUE) {
      for (i in 1:length(ellipse.pt)) {
        points(ellipse.pt[[i]], type="l", col=study.ellipse.col, lty=study.ellipse.lty,
               lwd=study.ellipse.lwd)
        }
    }
  
    points(x=ES[1], y=ES[2], pch=estimate.pch, cex=estimate.cex)
    if (estimate.ellipse.plot==TRUE) {
      points(ellipse(ACov, centre=ES), type="l", col=estimate.ellipse.col,
             lty=estimate.ellipse.lty, lwd=estimate.ellipse.lwd)
    }  
  # fixme
    if (randeff.ellipse.plot==TRUE) {
      points(ellipse(RE, centre=ES), type="l", col=randeff.ellipse.col,
             lty=randeff.ellipse.lty, lwd=randeff.ellipse.lwd)
    }
    if (univariate.plot==TRUE) {
      intervals <- c((1-interval)/2, (1+interval)/2)      
      x.m <- ES[1]
      x.se.limit <- qnorm(intervals, mean=x.m, sd=sqrt(ACov[1,1]))
      x.tau.limit <- qnorm(intervals, mean=x.m, sd=sqrt(RE[1,1]))
      abline(v=x.se.limit, col=univariate.lines.col, lty=univariate.lines.lty, lwd=univariate.lines.lwd)
      abline(v=x.tau.limit, col=univariate.lines.col, lty=univariate.lines.lty, lwd=univariate.lines.lwd)
      polygon(c(x.se.limit[1],x.m,x.se.limit[2],x.m),
              c(ylim[1], ylim[1]-0.02, ylim[1], ylim[1]+0.02), col=univariate.polygon.col)
      arrows(x.tau.limit[1], ylim[1], x.tau.limit[2], ylim[1], code=3, length=0.1,
             col=univariate.arrows.col, lwd=univariate.arrows.lwd)   
      y.m <- ES[2]
      y.se.limit <- qnorm(intervals, mean=y.m, sd=sqrt(ACov[2,2]))
      y.tau.limit <- qnorm(intervals, mean=y.m, sd=sqrt(RE[2,2]))
      abline(h=y.se.limit, col=univariate.lines.col, lty=univariate.lines.lty, lwd=univariate.lines.lwd)
      abline(h=y.tau.limit, col=univariate.lines.col, lty=univariate.lines.lty, lwd=univariate.lines.lwd)
      polygon(c(xlim[1], xlim[1]-0.02, xlim[1], xlim[1]+0.02),
              c(y.se.limit[1],y.m,y.se.limit[2],y.m),  col=univariate.polygon.col)
      arrows(xlim[1], y.tau.limit[1], xlim[1], y.tau.limit[2], code=3, length=0.1,
             col=univariate.arrows.col, lwd=univariate.arrows.lwd)
    }
  }
}


set.seed(1001)
x <- with(Berkey98, meta(y=cbind(PD, AL), v=cbind(var_PD, cov_PD_AL, var_AL)))
plot(x)
plot(x, effect.sizes=1:2)

my.df <- genNormMean(20, mu=c(0,0,0), sigma.between=vec2symMat(c(0.5,.2,.2,0.5,.2,0.5)),
                     sigma.within=vec2symMat(c(10,3,3,10,3,10)), n.mean=50, n.var=10, n.min=5)
my.meta <- meta(y=my.df[,4:6], v=my.df[,7:12])
meta1 <- rma(yi=my.df$y1, vi=my.df$v11)
meta2 <- rma(yi=my.df$y2, vi=my.df$v22)
meta3 <- rma(yi=my.df$y3, vi=my.df$v33)

plot(my.meta)
plot(my.meta, effect.sizes=c(2,3))
plot(my.meta, diag.panel=TRUE)
forest(meta1)
forest(meta2)
forest(meta3)
