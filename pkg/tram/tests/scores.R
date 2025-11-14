
library("mvtnorm")
library("tram")
library("numDeriv")

set.seed(2908)

N <- 100
J <- 4
Jp <- J * (J - 1) / 2
prm <- list(matrix(c(-1, 1, 1, -1), nrow = 2))
prm <- do.call("cbind", prm[rep(1, Jp)])[,1:Jp,drop = FALSE]

unif <- function(M, d)
    matrix(runif(M * d), ncol = d)

Genz <- function(M, d, fun = unif)
    list(w = t(fun(M, d = d)), M = M, fast = TRUE)

dgp <- function(N, J, conditional = FALSE, discrete = 0) {

    x <- runif(N)
    X <- cbind(1, x)
    lp <- X %*% prm
    L <- ltMatrices(t(lp), diag = FALSE, byrow = TRUE)
    if (!conditional) L <- invcholD(L)

    Z <- matrix(rnorm(J * N), ncol = N)
    Y <- solve(L, Z)
    d <- data.frame(t(Y), x = x)
    colnames(d)[1:J] <- paste0("Y", 1:J)
    if (length(discrete) == 1) discrete <- rep(discrete, J)
    if (any(discrete != 0)) 
        d[1:J] <- lapply(1:J, function(j) {
            if (discrete[j] == 0) return(d[[j]])
            dm <- abs(discrete[j])
            qn <- quantile(d[[j]], prob = c(1:(dm - 1)) / dm)
            ret <- cut(d[[j]], breaks = c(-Inf, qn, Inf), 
                       ordered_result = TRUE)
            if (discrete[j] > 0) return(ret)
            qn <- c(min(d[[j]]), qn, max(d[[j]]))
            ret <- (qn[-length(qn)] + diff(qn) / 2)[unclass(ret)]
            return(ret)
        })
    attr(d, "conditional") <- conditional
    attr(d, "discrete") <- discrete
    return(d)
}

### 
mlfun <- function(d, order = 1, formula = ~ x, args = NULL, probit = TRUE) {
    conditional <- attr(d, "conditional")
    discrete <- attr(d, "discrete")
    yvar <- colnames(d)[1:J]
    m <- lapply(yvar, function(y) {
        if (inherits(d[[y]], "factor")) {
            ret <- Polr(as.formula(paste(y, "~1")), data = d, method = "probit")
        } else {
            ret <- BoxCox(as.formula(paste(y, "~1")), data = d, order = order)
        }
        if (!probit) ret$model$todistr$name <- "CFG"
        ret
    })
    m$data <- d
    m$formula <- formula
    m$conditional <- conditional
    m$optim <- mltoptim(spg = list(maxit = 10000, 
                                   quiet = TRUE, ### checkGrad.tol = 1e-5,
                                   checkGrad = TRUE))["spg"]
    m$args <- args
    mm <- do.call("Mmlt", m)
    return(mm)
}

ORD <- 6
M <- 100

args <- NULL
if (J / 2 - 1 > 0)
    args <- Genz(M, J / 2 - 1)

d <- dgp(N = N, J = J, conditional = FALSE,
         discrete = rep(c(0, 10), each = J / 2))
mm <- mlfun(d, order = ORD, formula = ~ 1, args = args)
mm <- mlfun(d, order = ORD, formula = ~ x, args = args)

d <- dgp(N = N, J = J, conditional = TRUE,
         discrete = rep(c(0, 10), each = J / 2))
mm <- mlfun(d, order = ORD, formula = ~ 1, args = args)
mm <- mlfun(d, order = ORD, formula = ~ x, args = args)

