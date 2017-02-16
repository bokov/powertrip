dl<-function (pars, x)
{
    a<-pars[1];b<-pars[2];s<-pars[3];
    .expr2 <- exp(b * x)
    .expr3 <- a * .expr2
    .expr4 <- s * a
    .expr5 <- .expr2 - 1
    .expr6 <- .expr4 * .expr5
    .expr8 <- 1 + .expr6/b
    .expr9 <- s + 1
    .expr11 <- -.expr9/s
    .expr12 <- .expr8^.expr11
    .expr13 <- .expr3 * .expr12
    .expr16 <- .expr11 - 1
    .expr17 <- .expr8^.expr16
    .expr18 <- s * .expr5
    .expr19 <- .expr18/b
    .expr20 <- .expr11 * .expr19
    .expr21 <- .expr17 * .expr20
    .expr23 <- .expr2 * .expr12 + .expr3 * .expr21
    .expr25 <- .expr2 * .expr21
    .expr27 <- .expr8^(.expr16 - 1)
    .expr36 <- .expr13^2
    .expr39 <- .expr2 * x
    .expr41 <- .expr4 * .expr39
    .expr43 <- b^2
    .expr45 <- .expr41/b - .expr6/.expr43
    .expr46 <- .expr11 * .expr45
    .expr47 <- .expr17 * .expr46
    .expr50 <- a * .expr39
    .expr53 <- .expr27 * (.expr16 * .expr45)
    .expr68 <- .expr50 * .expr12 + .expr3 * .expr47
    .expr72 <- a * .expr5
    .expr73 <- .expr72/b
    .expr74 <- .expr11 * .expr73
    .expr76 <- log(.expr8)
    .expr78 <- s^2
    .expr80 <- 1/s - .expr9/.expr78
    .expr81 <- .expr76 * .expr80
    .expr83 <- .expr17 * .expr74 - .expr12 * .expr81
    .expr88 <- .expr27 * (.expr16 * .expr73) - .expr17 * .expr81
    .expr99 <- .expr3 * .expr83
    .expr104 <- .expr39 * x
    .expr107 <- .expr50 * .expr47
    .expr112 <- .expr41/.expr43
    .expr154 <- 1/.expr78
    .value <- log(.expr13)
    .grad <- array(0, c(length(.value), 3L), list(NULL, c("a",
        "b", "s")))
    .hessian <- array(0, c(length(.value), 3L, 3L), list(NULL,
        c("a", "b", "s"), c("a", "b", "s")))
    .grad[, "a"] <- .expr23/.expr13
    .hessian[, "a", "a"] <- (.expr25 + (.expr25 + .expr3 * (.expr27 *
        (.expr16 * .expr19) * .expr20)))/.expr13 - .expr23 *
        .expr23/.expr36
    .hessian[, "a", "b"] <- .hessian[, "b", "a"] <- (.expr39 *
        .expr12 + .expr2 * .expr47 + (.expr50 * .expr21 + .expr3 *
        (.expr53 * .expr20 + .expr17 * (.expr11 * (s * .expr39/b -
            .expr18/.expr43)))))/.expr13 - .expr23 * .expr68/.expr36
    .hessian[, "a", "s"] <- .hessian[, "s", "a"] <- (.expr2 *
        .expr83 + .expr3 * (.expr88 * .expr20 + .expr17 * (.expr11 *
        (.expr5/b) - .expr80 * .expr19)))/.expr13 - .expr23 *
        .expr99/.expr36
    .grad[, "b"] <- .expr68/.expr13
    .hessian[, "b", "b"] <- (a * .expr104 * .expr12 + .expr107 +
        (.expr107 + .expr3 * (.expr53 * .expr46 + .expr17 * (.expr11 *
            (.expr4 * .expr104/b - .expr112 - (.expr112 - .expr6 *
                (2 * b)/.expr43^2))))))/.expr13 - .expr68 * .expr68/.expr36
    .hessian[, "b", "s"] <- .hessian[, "s", "b"] <- (.expr50 *
        .expr83 + .expr3 * (.expr88 * .expr46 + .expr17 * (.expr11 *
        (.expr50/b - .expr72/.expr43) - .expr80 * .expr45)))/.expr13 -
        .expr68 * .expr99/.expr36
    .grad[, "s"] <- .expr99/.expr13
    .hessian[, "s", "s"] <- .expr3 * (.expr88 * .expr74 - .expr17 *
        (.expr80 * .expr73) - (.expr83 * .expr81 + .expr12 *
        (.expr73/.expr8 * .expr80 - .expr76 * (.expr154 + (.expr154 -
            .expr9 * (2 * s)/.expr78^2)))))/.expr13 - .expr99 *
        .expr99/.expr36
    attr(.value, "gradient") <- .grad
    attr(.value, "hessian") <- .hessian
    .value
}
