
library("mvtnorm")

### data generating process
dgp <- function(n = 100, 		### number of observations
                J = 2, 			### number of variables
                lambda = 2 * runif(J * (J - 1) / 2) - 1,
                                	### lower off-diagonal elements of Lambda
                df = 1, 		### chisq degrees of freedom for jth 
                                	### marginal distribution
                logchisq = FALSE,	### return(log(chisq))
                ncat = Inf		### bin observations
                )
{

    L <- ltMatrices(lambda, byrow = TRUE)
    Ls <- invcholD(L)
    Z <- matrix(rnorm(n * J), nrow = J)
    Zs <- solve(Ls, Z)
    df <- rep(df, length.out = J)
    Y <- do.call("data.frame", lapply(1:J, function(j) 
        qchisq(pnorm(Zs[j,], log = TRUE), df = df[j], log = TRUE)))

    logchisq <- rep(logchisq, length.out = J)
    if (any(logchisq)) Y[, logchisq] <- log(Y[,logchisq])

    if (any(is.finite(ncat))) {
        ncat <- rep(ncat, length.out = J)
        stopifnot(all(ncat >= 2) && isTRUE(all.equal(ncat, floor(ncat))))
        prob <- vector(mode = "list", length = J)
        prob[is.finite(ncat)] <- lapply(which(is.finite(ncat)), 
            function(j) sort(runif(ncat[j] - 1, min = .2, max = .8)))
        Y[,is.finite(ncat)] <- do.call("data.frame", 
            lapply(which(is.finite(ncat)), function(j) {
                qy <- quantile(Y[,j], prob = prob[[j]])
                yc <- cut(Y[,j], breaks = c(-Inf, qy, Inf), 
                          ordered_result = TRUE)
                ### remove empty levels
                yc <- yc[,drop = TRUE]
                if (nlevels(yc) == 2) 
                    class(yc) <- class(yc)[-1L]
                yc[, drop = TRUE]
            }))
    }

    colnames(Y) <- paste0("Y", 1:J)
    attr(Y, "L") <- L
    attr(Y, "positive") <- !logchisq
    return(Y)
}

summary(dgp(1000, ncat = 2))
summary(dgp(1000, ncat = 3))

summary(dgp(1000, J = 3, df = c(2, 1, 1), ncat = c(Inf, 4, 2)))

