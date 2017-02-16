dls<-function(pars,x)
{
    pars<-as.numeric(pars);
    a<-pars[1];b<-pars[2];s<-pars[3];
    .expr1 <- s * a
    .expr3 <- exp(b * x)
    .expr4 <- .expr3 - 1
    .expr5 <- .expr1 * .expr4
    .expr7 <- 1 + .expr5/b
    .expr9 <- -1/s
    .expr10 <- .expr7^.expr9
    .expr12 <- .expr9 - 1
    .expr13 <- .expr7^.expr12
    .expr14 <- s * .expr4
    .expr15 <- .expr14/b
    .expr16 <- .expr9 * .expr15
    .expr17 <- .expr13 * .expr16
    .expr20 <- .expr7^(.expr12 - 1)
    .expr26 <- .expr10^2
    .expr29 <- .expr3 * x
    .expr30 <- .expr1 * .expr29
    .expr32 <- b^2
    .expr34 <- .expr30/b - .expr5/.expr32
    .expr36 <- .expr20 * (.expr12 * .expr34)
    .expr46 <- .expr9 * .expr34
    .expr47 <- .expr13 * .expr46
    .expr51 <- a * .expr4
    .expr52 <- .expr51/b
    .expr55 <- log(.expr7)
    .expr56 <- s^2
    .expr57 <- 1/.expr56
    .expr58 <- .expr55 * .expr57
    .expr60 <- .expr20 * (.expr12 * .expr52) + .expr13 * .expr58
    .expr69 <- .expr9 * .expr52
    .expr72 <- .expr13 * .expr69 + .expr10 * .expr58
    .expr81 <- .expr30/.expr32
    .value <- log(.expr10)
    .grad <- array(0, c(length(.value), 3L), list(NULL, c("a",
        "b", "s")))
    .hessian <- array(0, c(length(.value), 3L, 3L), list(NULL,
        c("a", "b", "s"), c("a", "b", "s")))
    .grad[, "a"] <- .expr17/.expr10
    .hessian[, "a", "a"] <- .expr20 * (.expr12 * .expr15) * .expr16/.expr10 -
        .expr17 * .expr17/.expr26
    .hessian[, "a", "b"] <- .hessian[, "b", "a"] <- (.expr36 *
        .expr16 + .expr13 * (.expr9 * (s * .expr29/b - .expr14/.expr32)))/.expr10 -
        .expr17 * .expr47/.expr26
    .hessian[, "a", "s"] <- .hessian[, "s", "a"] <- (.expr60 *
        .expr16 + .expr13 * (.expr57 * .expr15 + .expr9 * (.expr4/b)))/.expr10 -
        .expr17 * .expr72/.expr26
    .grad[, "b"] <- .expr47/.expr10
    .hessian[, "b", "b"] <- (.expr36 * .expr46 + .expr13 * (.expr9 *
        (.expr1 * (.expr29 * x)/b - .expr81 - (.expr81 - .expr5 *
            (2 * b)/.expr32^2))))/.expr10 - .expr47 * .expr47/.expr26
    .hessian[, "b", "s"] <- .hessian[, "s", "b"] <- (.expr60 *
        .expr46 + .expr13 * (.expr57 * .expr34 + .expr9 * (a *
        .expr29/b - .expr51/.expr32)))/.expr10 - .expr47 * .expr72/.expr26
    .grad[, "s"] <- .expr72/.expr10
    .hessian[, "s", "s"] <- (.expr60 * .expr69 + .expr13 * (.expr57 *
        .expr52) + (.expr72 * .expr58 + .expr10 * (.expr52/.expr7 *
        .expr57 - .expr55 * (2 * s/.expr56^2))))/.expr10 - .expr72 *
        .expr72/.expr26
    attr(.value, "gradient") <- .grad
    attr(.value, "hessian") <- .hessian
    .value
}