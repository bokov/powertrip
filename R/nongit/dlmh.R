dlmh<-function(pars,x)
{
    pars<-as.numeric(pars);
    a<-pars[1];b<-pars[2];c<-pars[3];s<-pars[4];
    .expr2 <- exp(b * x)
    .expr3 <- a * .expr2
    .expr4 <- s * a
    .expr5 <- .expr2 - 1
    .expr6 <- .expr4 * .expr5
    .expr8 <- 1 + .expr6/b
    .expr10 <- c + .expr3/.expr8
    .expr13 <- s * .expr5
    .expr14 <- .expr13/b
    .expr15 <- .expr3 * .expr14
    .expr16 <- .expr8^2
    .expr18 <- .expr2/.expr8 - .expr15/.expr16
    .expr21 <- .expr2 * .expr14/.expr16
    .expr25 <- .expr16^2
    .expr31 <- .expr10^2
    .expr35 <- .expr2 * x
    .expr37 <- .expr4 * .expr35
    .expr39 <- b^2
    .expr41 <- .expr37/b - .expr6/.expr39
    .expr45 <- a * .expr35
    .expr55 <- 2 * (.expr41 * .expr8)
    .expr62 <- .expr3 * .expr41
    .expr64 <- .expr45/.expr8 - .expr62/.expr16
    .expr70 <- a * .expr5
    .expr71 <- .expr70/b
    .expr78 <- 2 * (.expr71 * .expr8)
    .expr84 <- .expr3 * .expr71
    .expr85 <- .expr84/.expr16
    .expr91 <- .expr35 * x
    .expr94 <- .expr45 * .expr41
    .expr99 <- .expr37/.expr39
    .value <- log(.expr10)
    .grad <- array(0, c(length(.value), 4L), list(NULL, c("a",
        "b", "c", "s")))
    .hessian <- array(0, c(length(.value), 4L, 4L), list(NULL,
        c("a", "b", "c", "s"), c("a", "b", "c", "s")))
    .grad[, "a"] <- .expr18/.expr10
    .hessian[, "a", "a"] <- -((.expr21 + (.expr21 - .expr15 *
        (2 * (.expr14 * .expr8))/.expr25))/.expr10 + .expr18 *
        .expr18/.expr31)
    .hessian[, "a", "b"] <- .hessian[, "b", "a"] <- (.expr35/.expr8 -
        .expr2 * .expr41/.expr16 - ((.expr45 * .expr14 + .expr3 *
        (s * .expr35/b - .expr13/.expr39))/.expr16 - .expr15 *
        .expr55/.expr25))/.expr10 - .expr18 * .expr64/.expr31
    .hessian[, "a", "c"] <- .hessian[, "c", "a"] <- -(.expr18/.expr31)
    .hessian[, "a", "s"] <- .hessian[, "s", "a"] <- -((.expr2 *
        .expr71/.expr16 + (.expr3 * (.expr5/b)/.expr16 - .expr15 *
        .expr78/.expr25))/.expr10 - .expr18 * .expr85/.expr31)
    .grad[, "b"] <- .expr64/.expr10
    .hessian[, "b", "b"] <- (a * .expr91/.expr8 - .expr94/.expr16 -
        ((.expr94 + .expr3 * (.expr4 * .expr91/b - .expr99 -
            (.expr99 - .expr6 * (2 * b)/.expr39^2)))/.expr16 -
            .expr62 * .expr55/.expr25))/.expr10 - .expr64 * .expr64/.expr31
    .hessian[, "b", "c"] <- .hessian[, "c", "b"] <- -(.expr64/.expr31)
    .hessian[, "b", "s"] <- .hessian[, "s", "b"] <- -((.expr45 *
        .expr71/.expr16 + (.expr3 * (.expr45/b - .expr70/.expr39)/.expr16 -
        .expr62 * .expr78/.expr25))/.expr10 - .expr64 * .expr85/.expr31)
    .grad[, "c"] <- 1/.expr10
    .hessian[, "c", "c"] <- -(1/.expr31)
    .hessian[, "c", "s"] <- .hessian[, "s", "c"] <- .expr85/.expr31
    .grad[, "s"] <- -(.expr85/.expr10)
    .hessian[, "s", "s"] <- .expr84 * .expr78/.expr25/.expr10 -
        .expr85 * .expr85/.expr31
    attr(.value, "gradient") <- .grad
    attr(.value, "hessian") <- .hessian
    .value
}