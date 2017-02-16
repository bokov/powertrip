dgh<-function(pars,x)
{
    pars<-as.numeric(pars);
    a<-pars[1];b<-pars[2];
    .expr2 <- exp(b * x)
    .expr3 <- a * .expr2
    .expr7 <- .expr3^2
    .expr10 <- .expr2 * x
    .expr12 <- a * .expr10
    .value <- log(.expr3)
    .grad <- array(0, c(length(.value), 2L), list(NULL, c("a",
        "b")))
    .hessian <- array(0, c(length(.value), 2L, 2L), list(NULL,
        c("a", "b"), c("a", "b")))
    .grad[, "a"] <- .expr2/.expr3
    .hessian[, "a", "a"] <- -(.expr2 * .expr2/.expr7)
    .hessian[, "a", "b"] <- .hessian[, "b", "a"] <- .expr10/.expr3 -
        .expr2 * .expr12/.expr7
    .grad[, "b"] <- .expr12/.expr3
    .hessian[, "b", "b"] <- a * (.expr10 * x)/.expr3 - .expr12 *
        .expr12/.expr7
    attr(.value, "gradient") <- .grad
    attr(.value, "hessian") <- .hessian
    .value
}