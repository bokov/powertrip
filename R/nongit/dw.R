dw<-function (pars, x) 
{
    pars<-as.numeric(pars);
    a<-pars[1];b<-pars[2];
    .expr1 <- a^b
    .expr2 <- .expr1 * b
    .expr3 <- b - 1
    .expr4 <- x^.expr3
    .expr5 <- .expr2 * .expr4
    .expr7 <- -a * x
    .expr8 <- .expr7^b
    .expr9 <- exp(.expr8)
    .expr10 <- .expr5 * .expr9
    .expr12 <- a^.expr3
    .expr13 <- .expr12 * b
    .expr14 <- .expr13 * b
    .expr15 <- .expr14 * .expr4
    .expr17 <- .expr7^.expr3
    .expr18 <- b * x
    .expr19 <- .expr17 * .expr18
    .expr20 <- .expr9 * .expr19
    .expr22 <- .expr15 * .expr9 - .expr5 * .expr20
    .expr24 <- .expr3 - 1
    .expr31 <- .expr15 * .expr20
    .expr45 <- .expr10^2
    .expr48 <- log(a)
    .expr55 <- log(x)
    .expr56 <- .expr4 * .expr55
    .expr60 <- log(.expr7)
    .expr61 <- .expr8 * .expr60
    .expr62 <- .expr9 * .expr61
    .expr65 <- .expr1 * .expr48
    .expr67 <- .expr65 * b + .expr1
    .expr70 <- .expr67 * .expr4 + .expr2 * .expr56
    .expr85 <- .expr70 * .expr9 + .expr5 * .expr62
    .expr95 <- .expr67 * .expr56
    .expr102 <- .expr70 * .expr62
    .value <- log(.expr10)
    .grad <- array(0, c(length(.value), 2L), list(NULL, c("a", 
        "b")))
    .hessian <- array(0, c(length(.value), 2L, 2L), list(NULL, 
        c("a", "b"), c("a", "b")))
    .grad[, "a"] <- .expr22/.expr10
    .hessian[, "a", "a"] <- (a^.expr24 * .expr3 * b * b * .expr4 * 
        .expr9 - .expr31 - (.expr31 - .expr5 * (.expr9 * (.expr7^.expr24 * 
        (.expr3 * x) * .expr18) + .expr20 * .expr19)))/.expr10 - 
        .expr22 * .expr22/.expr45
    .hessian[, "a", "b"] <- .hessian[, "b", "a"] <- ((((.expr12 * 
        .expr48 * b + .expr12) * b + .expr13) * .expr4 + .expr14 * 
        .expr56) * .expr9 + .expr15 * .expr62 - (.expr70 * .expr20 + 
        .expr5 * (.expr62 * .expr19 + .expr9 * (.expr17 * .expr60 * 
            .expr18 + .expr17 * x))))/.expr10 - .expr22 * .expr85/.expr45
    .grad[, "b"] <- .expr85/.expr10
    .hessian[, "b", "b"] <- (((.expr65 * .expr48 * b + .expr65 + 
        .expr65) * .expr4 + .expr95 + (.expr95 + .expr2 * (.expr56 * 
        .expr55))) * .expr9 + .expr102 + (.expr102 + .expr5 * 
        (.expr62 * .expr61 + .expr9 * (.expr61 * .expr60))))/.expr10 - 
        .expr85 * .expr85/.expr45
    attr(.value, "gradient") <- .grad
    attr(.value, "hessian") <- .hessian
    .value
}

