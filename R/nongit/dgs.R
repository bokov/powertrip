dgs<-function(pars,x)
{
    pars<-as.numeric(pars);
    a<-pars[1];b<-pars[2];
    .expr3 <- exp(b * x)
    .expr4 <- .expr3 - 1
    .expr5 <- -a * .expr4
    .expr9 <- .expr3 * x
    .expr11 <- b^2
    .expr15 <- a * .expr9
    .expr23 <- .expr15/.expr11
    .value <- .expr5/b
    .grad <- array(0, c(length(.value), 2L), list(NULL, c("a",
        "b")))
    .hessian <- array(0, c(length(.value), 2L, 2L), list(NULL,
        c("a", "b"), c("a", "b")))
    .grad[, "a"] <- -(.expr4/b)
    .hessian[, "a", "a"] <- 0
    .hessian[, "a", "b"] <- .hessian[, "b", "a"] <- -(.expr9/b -
        .expr4/.expr11)
    .grad[, "b"] <- -(.expr15/b + .expr5/.expr11)
    .hessian[, "b", "b"] <- -(a * (.expr9 * x)/b - .expr23 -
        (.expr23 + .expr5 * (2 * b)/.expr11^2))
    attr(.value, "gradient") <- .grad
    attr(.value, "hessian") <- .hessian
    .value
}