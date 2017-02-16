plot.summary.crqs<-function (x, nrow = 1, ncol = 1, CoxPHit = NULL, plcol='LightSkyBlue',ptcol='blue',
sigcol='red',lncol='blue',sigpch=20,sig=0.05,sigf=function(s){return(s)},...)
{
    args<-list(...)
    taus <- function(x) x$tau
    xx <- unlist(lapply(x, taus))
    coef <- lapply(x, coefficients)
    p <- nrow(coef[[1]])
    k <- ncol(coef[[1]])
    if (k != 6)
        stop("summary.crqs object has wrong column dimension")
    m <- length(xx)
    blab <- dimnames(coef[[1]])[[1]]
    a <- array(unlist(coef), c(p, k, m))
    if (length(CoxPHit))
        CoxQTE <- QTECox(CoxPHit)
    oldpar <- par(no.readonly = TRUE)
    par(mfrow = c(nrow, ncol))
    for (i in 2:p) {
        b <- a[i, 1, ]
        bl <- a[i, 2, ]
        bu <- a[i, 3, ]
	s <- a[i, 6, ]
	sadj <- sigf(s);
	ylim <- range(bl,bu,args$ylim)
        plot(rep(xx, 2), c(bl, bu), xlab = "", ylab = "", type = "n", ylim = ylim)
        title(ifelse(is.null(args$main),paste(blab[i]),args$main), cex = ifelse(is.null(args$cex),.75,args$cex))
        polygon(c(xx, rev(xx)), c(bl, rev(bu)), col = plcol)
        points(xx, b, cex = ifelse(is.null(args$cex),0.5,0.67*args$cex), pch = ifelse(is.null(args$pch),'o',args$pch), col = ptcol)
        lines(xx, b, col = lncol)
	points(xx[sadj<sig], b[sadj<sig], cex = ifelse(is.null(args$cex),0.6,0.8*args$cex), lwd=1, pch = sigpch, col = sigcol)
        abline(h = 0)
        if (length(CoxPHit)) {
            lines(CoxQTE$taus, CoxQTE$QTE[, i - 1], col = "red")
        }
    }
    par(oldpar)
}
