dg<-function (pars, x)
{
    a<-pars[1];b<-pars[2];
    .expr1 <- b * x
    .expr2 <- exp(.expr1)
    .expr3 <- .expr2 - 1
    .expr4 <- a * .expr3
    .expr7 <- exp(.expr1 - .expr4/b)
    .expr8 <- a * .expr7
    .expr10 <- .expr3/b
    .expr11 <- .expr7 * .expr10
    .expr13 <- .expr7 - a * .expr11
    .expr21 <- .expr8^2
    .expr25 <- .expr2 * x
    .expr26 <- a * .expr25
    .expr28 <- b^2
    .expr31 <- x - (.expr26/b - .expr4/.expr28)
    .expr32 <- .expr7 * .expr31
    .expr42 <- a * .expr32
    .expr51 <- .expr26/.expr28
    .value <- log(.expr8)
    .grad <- array(0, c(length(.value), 2L), list(NULL, c("a",
        "b")))
    .hessian <- array(0, c(length(.value), 2L, 2L), list(NULL,
        c("a", "b"), c("a", "b")))
    .grad[, "a"] <- .expr13/.expr8
    .hessian[, "a", "a"] <- -((.expr11 + (.expr11 - a * (.expr11 *
        .expr10)))/.expr8 + .expr13 * .expr13/.expr21)
    .hessian[, "a", "b"] <- .hessian[, "b", "a"] <- (.expr32 -
        a * (.expr32 * .expr10 + .expr7 * (.expr25/b - .expr3/.expr28)))/.expr8 -
        .expr13 * .expr42/.expr21
    .grad[, "b"] <- .expr42/.expr8
    .hessian[, "b", "b"] <- a * (.expr32 * .expr31 - .expr7 *
        (a * (.expr25 * x)/b - .expr51 - (.expr51 - .expr4 *
            (2 * b)/.expr28^2)))/.expr8 - .expr42 * .expr42/.expr21
    attr(.value, "gradient") <- .grad
    attr(.value, "hessian") <- .hessian
    .value
}
