dlms<-function(pars,x)
{
    pars<-as.numeric(pars);
    a<-pars[1];b<-pars[2];c<-pars[3];s<-pars[4];
    .expr3 <- exp(-c * x)
    .expr4 <- s * a
    .expr6 <- exp(b * x)
    .expr7 <- .expr6 - 1
    .expr8 <- .expr4 * .expr7
    .expr10 <- 1 + .expr8/b
    .expr12 <- -1/s
    .expr13 <- .expr10^.expr12
    .expr14 <- .expr3 * .expr13
    .expr16 <- .expr12 - 1
    .expr17 <- .expr10^.expr16
    .expr18 <- s * .expr7
    .expr19 <- .expr18/b
    .expr20 <- .expr12 * .expr19
    .expr21 <- .expr17 * .expr20
    .expr22 <- .expr3 * .expr21
    .expr25 <- .expr10^(.expr16 - 1)
    .expr32 <- .expr14^2
    .expr35 <- .expr6 * x
    .expr36 <- .expr4 * .expr35
    .expr38 <- b^2
    .expr40 <- .expr36/b - .expr8/.expr38
    .expr42 <- .expr25 * (.expr16 * .expr40)
    .expr53 <- .expr12 * .expr40
    .expr54 <- .expr17 * .expr53
    .expr55 <- .expr3 * .expr54
    .expr59 <- .expr3 * x
    .expr62 <- .expr59 * .expr13
    .expr67 <- a * .expr7
    .expr68 <- .expr67/b
    .expr71 <- log(.expr10)
    .expr72 <- s^2
    .expr73 <- 1/.expr72
    .expr74 <- .expr71 * .expr73
    .expr76 <- .expr25 * (.expr16 * .expr68) + .expr17 * .expr74
    .expr86 <- .expr12 * .expr68
    .expr89 <- .expr17 * .expr86 + .expr13 * .expr74
    .expr90 <- .expr3 * .expr89
    .expr99 <- .expr36/.expr38
    .value <- log(.expr14)
    .grad <- array(0, c(length(.value), 4L), list(NULL, c("a",
        "b", "c", "s")))
    .hessian <- array(0, c(length(.value), 4L, 4L), list(NULL,
        c("a", "b", "c", "s"), c("a", "b", "c", "s")))
    .grad[, "a"] <- .expr22/.expr14
    .hessian[, "a", "a"] <- .expr3 * (.expr25 * (.expr16 * .expr19) *
        .expr20)/.expr14 - .expr22 * .expr22/.expr32
    .hessian[, "a", "b"] <- .hessian[, "b", "a"] <- .expr3 *
        (.expr42 * .expr20 + .expr17 * (.expr12 * (s * .expr35/b -
            .expr18/.expr38)))/.expr14 - .expr22 * .expr55/.expr32
    .hessian[, "a", "c"] <- .hessian[, "c", "a"] <- -(.expr59 *
        .expr21/.expr14 - .expr22 * .expr62/.expr32)
    .hessian[, "a", "s"] <- .hessian[, "s", "a"] <- .expr3 *
        (.expr76 * .expr20 + .expr17 * (.expr73 * .expr19 + .expr12 *
            (.expr7/b)))/.expr14 - .expr22 * .expr90/.expr32
    .grad[, "b"] <- .expr55/.expr14
    .hessian[, "b", "b"] <- .expr3 * (.expr42 * .expr53 + .expr17 *
        (.expr12 * (.expr4 * (.expr35 * x)/b - .expr99 - (.expr99 -
            .expr8 * (2 * b)/.expr38^2))))/.expr14 - .expr55 *
        .expr55/.expr32
    .hessian[, "b", "c"] <- .hessian[, "c", "b"] <- -(.expr59 *
        .expr54/.expr14 - .expr55 * .expr62/.expr32)
    .hessian[, "b", "s"] <- .hessian[, "s", "b"] <- .expr3 *
        (.expr76 * .expr53 + .expr17 * (.expr73 * .expr40 + .expr12 *
            (a * .expr35/b - .expr67/.expr38)))/.expr14 - .expr55 *
        .expr90/.expr32
    .grad[, "c"] <- -(.expr62/.expr14)
    .hessian[, "c", "c"] <- .expr59 * x * .expr13/.expr14 - .expr62 *
        .expr62/.expr32
    .hessian[, "c", "s"] <- .hessian[, "s", "c"] <- -(.expr59 *
        .expr89/.expr14 - .expr62 * .expr90/.expr32)
    .grad[, "s"] <- .expr90/.expr14
    .hessian[, "s", "s"] <- .expr3 * (.expr76 * .expr86 + .expr17 *
        (.expr73 * .expr68) + (.expr89 * .expr74 + .expr13 *
        (.expr68/.expr10 * .expr73 - .expr71 * (2 * s/.expr72^2))))/.expr14 -
        .expr90 * .expr90/.expr32
    attr(.value, "gradient") <- .grad
    attr(.value, "hessian") <- .hessian
    .value
}