dwh<-function (pars, x)
{
    pars<-as.numeric(pars);
    a<-pars[1];b<-pars[2];
    .expr1 <- a^b
    .expr2 <- .expr1 * b
    .expr3 <- b - 1
    .expr4 <- x^.expr3
    .expr5 <- .expr2 * .expr4
    .expr7 <- a^.expr3
    .expr8 <- .expr7 * b
    .expr9 <- .expr8 * b
    .expr10 <- .expr9 * .expr4
    .expr20 <- .expr5^2
    .expr23 <- log(a)
    .expr30 <- log(x)
    .expr31 <- .expr4 * .expr30
    .expr35 <- .expr1 * .expr23
    .expr37 <- .expr35 * b + .expr1
    .expr40 <- .expr37 * .expr4 + .expr2 * .expr31
    .expr50 <- .expr37 * .expr31
    .value <- log(.expr5)
    .grad <- array(0, c(length(.value), 2L), list(NULL, c("a",
        "b")))
    .hessian <- array(0, c(length(.value), 2L, 2L), list(NULL,
        c("a", "b"), c("a", "b")))
    .grad[, "a"] <- .expr10/.expr5
    .hessian[, "a", "a"] <- a^(.expr3 - 1) * .expr3 * b * b *
        .expr4/.expr5 - .expr10 * .expr10/.expr20
    .hessian[, "a", "b"] <- .hessian[, "b", "a"] <- (((.expr7 *
        .expr23 * b + .expr7) * b + .expr8) * .expr4 + .expr9 *
        .expr31)/.expr5 - .expr10 * .expr40/.expr20
    .grad[, "b"] <- .expr40/.expr5
    .hessian[, "b", "b"] <- ((.expr35 * .expr23 * b + .expr35 +
        .expr35) * .expr4 + .expr50 + (.expr50 + .expr2 * (.expr31 *
        .expr30)))/.expr5 - .expr40 * .expr40/.expr20
    attr(.value, "gradient") <- .grad
    attr(.value, "hessian") <- .hessian
    .value
}
