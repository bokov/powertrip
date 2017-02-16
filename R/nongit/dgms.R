dgms<-function(pars,x)
{
    pars<-as.numeric(pars);
    a<-pars[1];b<-pars[2];c<-pars[3];
    .expr3 <- exp(b * x)
    .expr4 <- .expr3 - 1
    .expr5 <- -a * .expr4
    .expr11 <- .expr3 * x
    .expr13 <- b^2
    .expr17 <- a * .expr11
    .expr25 <- .expr17/.expr13
    .value <- .expr5/b - c * x
    .grad <- array(0, c(length(.value), 3L), list(NULL, c("a",
        "b", "c")))
    .hessian <- array(0, c(length(.value), 3L, 3L), list(NULL,
        c("a", "b", "c"), c("a", "b", "c")))
    .grad[, "a"] <- -(.expr4/b)
    .hessian[, "a", "a"] <- 0
    .hessian[, "a", "b"] <- .hessian[, "b", "a"] <- -(.expr11/b -
        .expr4/.expr13)
    .hessian[, "a", "c"] <- .hessian[, "c", "a"] <- 0
    .grad[, "b"] <- -(.expr17/b + .expr5/.expr13)
    .hessian[, "b", "b"] <- -(a * (.expr11 * x)/b - .expr25 -
        (.expr25 + .expr5 * (2 * b)/.expr13^2))
    .hessian[, "b", "c"] <- .hessian[, "c", "b"] <- 0
    .grad[, "c"] <- -x
    .hessian[, "c", "c"] <- 0
    attr(.value, "gradient") <- .grad
    attr(.value, "hessian") <- .hessian
    .value
}