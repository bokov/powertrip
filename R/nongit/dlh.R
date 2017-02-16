dlh<-function(pars,x)
{
    pars<-as.numeric(pars);
    a<-pars[1];b<-pars[2];s<-pars[3];
    .expr2 <- exp(b * x)
    .expr3 <- a * .expr2
    .expr4 <- s * a
    .expr5 <- .expr2 - 1
    .expr6 <- .expr4 * .expr5
    .expr8 <- 1 + .expr6/b
    .expr9 <- .expr3/.expr8
    .expr12 <- s * .expr5
    .expr13 <- .expr12/b
    .expr14 <- .expr3 * .expr13
    .expr15 <- .expr8^2
    .expr17 <- .expr2/.expr8 - .expr14/.expr15
    .expr20 <- .expr2 * .expr13/.expr15
    .expr24 <- .expr15^2
    .expr30 <- .expr9^2
    .expr34 <- .expr2 * x
    .expr36 <- .expr4 * .expr34
    .expr38 <- b^2
    .expr40 <- .expr36/b - .expr6/.expr38
    .expr44 <- a * .expr34
    .expr54 <- 2 * (.expr40 * .expr8)
    .expr61 <- .expr3 * .expr40
    .expr63 <- .expr44/.expr8 - .expr61/.expr15
    .expr67 <- a * .expr5
    .expr68 <- .expr67/b
    .expr75 <- 2 * (.expr68 * .expr8)
    .expr81 <- .expr3 * .expr68
    .expr82 <- .expr81/.expr15
    .expr88 <- .expr34 * x
    .expr91 <- .expr44 * .expr40
    .expr96 <- .expr36/.expr38
    .value <- log(.expr9)
    .grad <- array(0, c(length(.value), 3L), list(NULL, c("a",
        "b", "s")))
    .hessian <- array(0, c(length(.value), 3L, 3L), list(NULL,
        c("a", "b", "s"), c("a", "b", "s")))
    .grad[, "a"] <- .expr17/.expr9
    .hessian[, "a", "a"] <- -((.expr20 + (.expr20 - .expr14 *
        (2 * (.expr13 * .expr8))/.expr24))/.expr9 + .expr17 *
        .expr17/.expr30)
    .hessian[, "a", "b"] <- .hessian[, "b", "a"] <- (.expr34/.expr8 -
        .expr2 * .expr40/.expr15 - ((.expr44 * .expr13 + .expr3 *
        (s * .expr34/b - .expr12/.expr38))/.expr15 - .expr14 *
        .expr54/.expr24))/.expr9 - .expr17 * .expr63/.expr30
    .hessian[, "a", "s"] <- .hessian[, "s", "a"] <- -((.expr2 *
        .expr68/.expr15 + (.expr3 * (.expr5/b)/.expr15 - .expr14 *
        .expr75/.expr24))/.expr9 - .expr17 * .expr82/.expr30)
    .grad[, "b"] <- .expr63/.expr9
    .hessian[, "b", "b"] <- (a * .expr88/.expr8 - .expr91/.expr15 -
        ((.expr91 + .expr3 * (.expr4 * .expr88/b - .expr96 -
            (.expr96 - .expr6 * (2 * b)/.expr38^2)))/.expr15 -
            .expr61 * .expr54/.expr24))/.expr9 - .expr63 * .expr63/.expr30
    .hessian[, "b", "s"] <- .hessian[, "s", "b"] <- -((.expr44 *
        .expr68/.expr15 + (.expr3 * (.expr44/b - .expr67/.expr38)/.expr15 -
        .expr61 * .expr75/.expr24))/.expr9 - .expr63 * .expr82/.expr30)
    .grad[, "s"] <- -(.expr82/.expr9)
    .hessian[, "s", "s"] <- .expr81 * .expr75/.expr24/.expr9 -
        .expr82 * .expr82/.expr30
    attr(.value, "gradient") <- .grad
    attr(.value, "hessian") <- .hessian
    .value
}