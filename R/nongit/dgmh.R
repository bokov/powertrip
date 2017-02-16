dgmh<-function(pars,x)
{
    pars<-as.numeric(pars);
    a<-pars[1];b<-pars[2];c<-pars[3];
    .expr2 <- exp(b * x)
    .expr4 <- c + a * .expr2
    .expr8 <- .expr4^2
    .expr11 <- .expr2 * x
    .expr13 <- a * .expr11
    .value <- log(.expr4)
    .grad <- array(0, c(length(.value), 3L), list(NULL, c("a",
        "b", "c")))
    .hessian <- array(0, c(length(.value), 3L, 3L), list(NULL,
        c("a", "b", "c"), c("a", "b", "c")))
    .grad[, "a"] <- .expr2/.expr4
    .hessian[, "a", "a"] <- -(.expr2 * .expr2/.expr8)
    .hessian[, "a", "b"] <- .hessian[, "b", "a"] <- .expr11/.expr4 -
        .expr2 * .expr13/.expr8
    .hessian[, "a", "c"] <- .hessian[, "c", "a"] <- -(.expr2/.expr8)
    .grad[, "b"] <- .expr13/.expr4
    .hessian[, "b", "b"] <- a * (.expr11 * x)/.expr4 - .expr13 *
        .expr13/.expr8
    .hessian[, "b", "c"] <- .hessian[, "c", "b"] <- -(.expr13/.expr8)
    .grad[, "c"] <- 1/.expr4
    .hessian[, "c", "c"] <- -(1/.expr8)
    attr(.value, "gradient") <- .grad
    attr(.value, "hessian") <- .hessian
    .value
}