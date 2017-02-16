library(survival);

doquant <- function (p, time, surv, upper, lower) 
{
    findq <- function(x, y, p) {
        if (all(is.na(y))) 
            return(rep(NA, length(p)))
        toss <- duplicated(y)
        if (any(toss)) {
            newy <- y[!toss]
            newx <- x[!toss]
        }
        else {
            newy <- y
            newx <- x
        }
        indx <- approx(c(0, 1 - newy), c(0:length(newy)), p)$y
        indx2 <- ceiling(indx)
        result <- newx[indx2]
        if (any(!is.na(indx) & indx == indx2)) {
            special <- which(indx == indx2)
            upper <- c(newx, max(x))[indx2[special] + 1]
            result[special] <- (result[special] + upper)/2
        }
        result
    }
    qq <- findq(time, surv, p)
    if (missing(upper)) 
        qq
    else rbind(qq, findq(time, lower, p), findq(time, upper, 
        p))
}

qsf <- survival:::quantile.survfit;
environment(qsf) <- .GlobalEnv;
