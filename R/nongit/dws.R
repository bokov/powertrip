dws<-function(pars,x)
{
    pars<-as.numeric(pars);
    a<-pars[1];b<-pars[2];
    .expr1 <- a * x
    .expr2 <- .expr1^b
    .expr4 <- b - 1
    .expr5 <- .expr1^.expr4
    .expr6 <- b * x
    .expr15 <- log(.expr1)
    .expr21 <- .expr2 * .expr15
    .value <- -.expr2
    .grad <- array(0, c(length(.value), 2L), list(NULL, c("a",
        "b")))
    .hessian <- array(0, c(length(.value), 2L, 2L), list(NULL,
        c("a", "b"), c("a", "b")))
    .grad[, "a"] <- -(.expr5 * .expr6)
    .hessian[, "a", "a"] <- -(.expr1^(.expr4 - 1) * (.expr4 *
        x) * .expr6)
    .hessian[, "a", "b"] <- .hessian[, "b", "a"] <- -(.expr5 *
        .expr15 * .expr6 + .expr5 * x)
    .grad[, "b"] <- -.expr21
    .hessian[, "b", "b"] <- -(.expr21 * .expr15)
    attr(.value, "gradient") <- .grad
    attr(.value, "hessian") <- .hessian
    .value
}