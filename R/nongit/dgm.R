dgm<-function (pars, x)
{
    a<-pars[1];b<-pars[2];c<-pars[3];
    .expr2 <- exp(b * x)
    .expr4 <- c + a * .expr2
    .expr6 <- .expr2 - 1
    .expr7 <- -a * .expr6
    .expr11 <- exp(.expr7/b - c * x)
    .expr12 <- .expr4 * .expr11
    .expr15 <- .expr6/b
    .expr16 <- .expr11 * .expr15
    .expr18 <- .expr2 * .expr11 - .expr4 * .expr16
    .expr20 <- .expr2 * .expr16
    .expr27 <- .expr12^2
    .expr31 <- .expr2 * x
    .expr33 <- a * .expr31
    .expr35 <- b^2
    .expr37 <- .expr33/b + .expr7/.expr35
    .expr38 <- .expr11 * .expr37
    .expr54 <- .expr33 * .expr11 - .expr4 * .expr38
    .expr58 <- .expr11 * x
    .expr66 <- .expr11 - .expr4 * .expr58
    .expr73 <- a * (.expr31 * x)
    .expr75 <- .expr33 * .expr38
    .expr78 <- .expr33/.expr35
    .value <- log(.expr12)
    .grad <- array(0, c(length(.value), 3L), list(NULL, c("a",
        "b", "c")))
    .hessian <- array(0, c(length(.value), 3L, 3L), list(NULL,
        c("a", "b", "c"), c("a", "b", "c")))
    .grad[, "a"] <- .expr18/.expr12
    .hessian[, "a", "a"] <- -((.expr20 + (.expr20 - .expr4 *
        (.expr16 * .expr15)))/.expr12 + .expr18 * .expr18/.expr27)
    .hessian[, "a", "b"] <- .hessian[, "b", "a"] <- (.expr31 *
        .expr11 - .expr2 * .expr38 - (.expr33 * .expr16 + .expr4 *
        (.expr11 * (.expr31/b - .expr6/.expr35) - .expr38 * .expr15)))/.expr12 -
        .expr18 * .expr54/.expr27
    .hessian[, "a", "c"] <- .hessian[, "c", "a"] <- -((.expr2 *
        .expr58 + (.expr16 - .expr4 * (.expr58 * .expr15)))/.expr12 +
        .expr18 * .expr66/.expr27)
    .grad[, "b"] <- .expr54/.expr12
    .hessian[, "b", "b"] <- (.expr73 * .expr11 - .expr75 - (.expr75 +
        .expr4 * (.expr11 * (.expr73/b - .expr78 - (.expr78 +
            .expr7 * (2 * b)/.expr35^2)) - .expr38 * .expr37)))/.expr12 -
        .expr54 * .expr54/.expr27
    .hessian[, "b", "c"] <- .hessian[, "c", "b"] <- -((.expr33 *
        .expr58 + (.expr38 - .expr4 * (.expr58 * .expr37)))/.expr12 +
        .expr54 * .expr66/.expr27)
    .grad[, "c"] <- .expr66/.expr12
    .hessian[, "c", "c"] <- -((.expr58 + (.expr58 - .expr4 *
        (.expr58 * x)))/.expr12 + .expr66 * .expr66/.expr27)
    attr(.value, "gradient") <- .grad
    attr(.value, "hessian") <- .hessian
    .value
}
