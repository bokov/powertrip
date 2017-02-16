dlm<-function (pars, x)
{
    a<-pars[1];b<-pars[2];c<-pars[3];s<-pars[4];
    .expr2 <- exp(b * x)
    .expr3 <- a * .expr2
    .expr4 <- s * a
    .expr5 <- .expr2 - 1
    .expr6 <- .expr4 * .expr5
    .expr8 <- 1 + .expr6/b
    .expr10 <- c + .expr3/.expr8
    .expr13 <- exp(-c * x)
    .expr15 <- -1/s
    .expr16 <- .expr8^.expr15
    .expr17 <- .expr13 * .expr16
    .expr18 <- .expr10 * .expr17
    .expr21 <- s * .expr5
    .expr22 <- .expr21/b
    .expr23 <- .expr3 * .expr22
    .expr24 <- .expr8^2
    .expr26 <- .expr2/.expr8 - .expr23/.expr24
    .expr28 <- .expr15 - 1
    .expr29 <- .expr8^.expr28
    .expr30 <- .expr15 * .expr22
    .expr31 <- .expr29 * .expr30
    .expr32 <- .expr13 * .expr31
    .expr34 <- .expr26 * .expr17 + .expr10 * .expr32
    .expr36 <- .expr26 * .expr32
    .expr38 <- .expr2 * .expr22/.expr24
    .expr42 <- .expr24^2
    .expr49 <- .expr8^(.expr28 - 1)
    .expr59 <- .expr18^2
    .expr62 <- .expr2 * x
    .expr64 <- .expr4 * .expr62
    .expr66 <- b^2
    .expr68 <- .expr64/b - .expr6/.expr66
    .expr72 <- a * .expr62
    .expr77 <- s * .expr62/b - .expr21/.expr66
    .expr82 <- 2 * (.expr68 * .expr8)
    .expr88 <- .expr15 * .expr68
    .expr89 <- .expr29 * .expr88
    .expr90 <- .expr13 * .expr89
    .expr94 <- .expr3 * .expr68
    .expr96 <- .expr72/.expr8 - .expr94/.expr24
    .expr99 <- .expr49 * (.expr28 * .expr68)
    .expr111 <- .expr96 * .expr17 + .expr10 * .expr90
    .expr115 <- .expr13 * x
    .expr119 <- .expr115 * .expr16
    .expr124 <- .expr17 - .expr10 * .expr119
    .expr128 <- a * .expr5
    .expr129 <- .expr128/b
    .expr130 <- .expr15 * .expr129
    .expr132 <- log(.expr8)
    .expr133 <- s^2
    .expr134 <- 1/.expr133
    .expr135 <- .expr132 * .expr134
    .expr137 <- .expr29 * .expr130 + .expr16 * .expr135
    .expr138 <- .expr13 * .expr137
    .expr142 <- .expr5/b
    .expr146 <- 2 * (.expr129 * .expr8)
    .expr156 <- .expr49 * (.expr28 * .expr129) + .expr29 * .expr135
    .expr165 <- .expr3 * .expr129
    .expr166 <- .expr165/.expr24
    .expr173 <- .expr10 * .expr138 - .expr166 * .expr17
    .expr178 <- .expr62 * x
    .expr181 <- .expr72 * .expr68
    .expr186 <- .expr64/.expr66
    .expr193 <- .expr4 * .expr178/b - .expr186 - (.expr186 -
        .expr6 * (2 * b)/.expr66^2)
    .expr202 <- .expr96 * .expr90
    .expr230 <- .expr72/b - .expr128/.expr66
    .expr292 <- .expr166 * .expr138
    .value <- log(.expr18)
    .grad <- array(0, c(length(.value), 4L), list(NULL, c("a",
        "b", "c", "s")))
    .hessian <- array(0, c(length(.value), 4L, 4L), list(NULL,
        c("a", "b", "c", "s"), c("a", "b", "c", "s")))
    .grad[, "a"] <- .expr34/.expr18
    .hessian[, "a", "a"] <- (.expr36 - (.expr38 + (.expr38 -
        .expr23 * (2 * (.expr22 * .expr8))/.expr42)) * .expr17 +
        (.expr36 + .expr10 * (.expr13 * (.expr49 * (.expr28 *
            .expr22) * .expr30))))/.expr18 - .expr34 * .expr34/.expr59
    .hessian[, "a", "b"] <- .hessian[, "b", "a"] <- ((.expr62/.expr8 -
        .expr2 * .expr68/.expr24 - ((.expr72 * .expr22 + .expr3 *
        .expr77)/.expr24 - .expr23 * .expr82/.expr42)) * .expr17 +
        .expr26 * .expr90 + (.expr96 * .expr32 + .expr10 * (.expr13 *
        (.expr99 * .expr30 + .expr29 * (.expr15 * .expr77)))))/.expr18 -
        .expr34 * .expr111/.expr59
    .hessian[, "a", "c"] <- .hessian[, "c", "a"] <- (.expr32 -
        .expr10 * (.expr115 * .expr31) - .expr26 * .expr119)/.expr18 -
        .expr34 * .expr124/.expr59
    .hessian[, "a", "s"] <- .hessian[, "s", "a"] <- (.expr26 *
        .expr138 - (.expr2 * .expr129/.expr24 + (.expr3 * .expr142/.expr24 -
        .expr23 * .expr146/.expr42)) * .expr17 + (.expr10 * (.expr13 *
        (.expr156 * .expr30 + .expr29 * (.expr134 * .expr22 +
            .expr15 * .expr142))) - .expr166 * .expr32))/.expr18 -
        .expr34 * .expr173/.expr59
    .grad[, "b"] <- .expr111/.expr18
    .hessian[, "b", "b"] <- ((a * .expr178/.expr8 - .expr181/.expr24 -
        ((.expr181 + .expr3 * .expr193)/.expr24 - .expr94 * .expr82/.expr42)) *
        .expr17 + .expr202 + (.expr202 + .expr10 * (.expr13 *
        (.expr99 * .expr88 + .expr29 * (.expr15 * .expr193)))))/.expr18 -
        .expr111 * .expr111/.expr59
    .hessian[, "b", "c"] <- .hessian[, "c", "b"] <- (.expr90 -
        .expr10 * (.expr115 * .expr89) - .expr96 * .expr119)/.expr18 -
        .expr111 * .expr124/.expr59
    .hessian[, "b", "s"] <- .hessian[, "s", "b"] <- (.expr96 *
        .expr138 - (.expr72 * .expr129/.expr24 + (.expr3 * .expr230/.expr24 -
        .expr94 * .expr146/.expr42)) * .expr17 + (.expr10 * (.expr13 *
        (.expr156 * .expr88 + .expr29 * (.expr134 * .expr68 +
            .expr15 * .expr230))) - .expr166 * .expr90))/.expr18 -
        .expr111 * .expr173/.expr59
    .grad[, "c"] <- .expr124/.expr18
    .hessian[, "c", "c"] <- -((.expr119 + (.expr119 - .expr10 *
        (.expr115 * x * .expr16)))/.expr18 + .expr124 * .expr124/.expr59)
    .hessian[, "c", "s"] <- .hessian[, "s", "c"] <- (.expr138 -
        (.expr10 * (.expr115 * .expr137) - .expr166 * .expr119))/.expr18 -
        .expr124 * .expr173/.expr59
    .grad[, "s"] <- .expr173/.expr18
    .hessian[, "s", "s"] <- (.expr10 * (.expr13 * (.expr156 *
        .expr130 + .expr29 * (.expr134 * .expr129) + (.expr137 *
        .expr135 + .expr16 * (.expr129/.expr8 * .expr134 - .expr132 *
        (2 * s/.expr133^2))))) - .expr292 - (.expr292 - .expr165 *
        .expr146/.expr42 * .expr17))/.expr18 - .expr173 * .expr173/.expr59
    attr(.value, "gradient") <- .grad
    attr(.value, "hessian") <- .hessian
    .value
}
