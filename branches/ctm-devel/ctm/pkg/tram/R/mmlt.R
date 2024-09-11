
Mmlt <- function(..., formula = ~ 1, data, conditional = FALSE, 
                 theta = NULL, fixed = NULL, scale = FALSE,
                 optim = mltoptim(auglag = list(maxtry = 5)), args = list(seed = 1, M = 1000), 
                 dofit = TRUE, domargins = TRUE, sequentialfit = FALSE)
{

    call <- match.call()  
    ret <- mmlt(..., formula = formula, data = data, conditional = conditional, 
                theta = theta, fixed = fixed, scale = scale,
                optim = optim, args = args,
                dofit = dofit, domargins = domargins)
                ###, sequentialfit = FALSE)
    ret$call <- call
    class(ret) <- c("Mmlt", class(ret))
    ret

}

summary.Mmlt <- function(object, ...) {
    ret <- list(call = object$call,
                #                tram = object$tram,
                test = cftest(object, 
                              parm = names(coef(object, with_baseline = FALSE))),
                ll = logLik(object))
    class(ret) <- "summary.Mmlt"
    ret
}

print.summary.Mmlt <- function(x, digits = max(3L, getOption("digits") - 3L), 
                               ...) {
    cat("\n", x$mmlt, "\n")
    cat("\nCall:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    pq <- x$test$test
    mtests <- cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
    colnames(mtests) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    sig <- .Machine$double.eps
    printCoefmat(mtests, digits = digits, has.Pvalue = TRUE, 
                 P.values = TRUE, eps.Pvalue = sig)
    cat("\nLog-Likelihood:\n ", x$ll, " (df = ", attr(x$ll, "df"), 
          ")", sep = "")
    cat("\n\n")
    invisible(x)
}

    
confregion <- function(object, level = .95, ...)
    UseMethod("confregion")

confregion.Mmlt <- function(object, level = .95, newdata, K = 250, ...) {

    if (!missing(newdata)) stopifnot(nrow(newdata) == 1)

    Linv <- coef(object, newdata = newdata, type = "Lambdainv")
    Linv <- as.array(Linv)[,,1]
    J <- nrow(Linv)

    q <- qchisq(level, df = J)

    if (J == 2) {
        angle <- seq(0, 2 * pi, length = K)
        x <- cbind(cos(angle), sin(angle))
    } else {
        x <- matrix(rnorm(K * J), nrow = K, ncol = J)
        x <- x / sqrt(rowSums(x^2))
    }
    x <- sqrt(q) * x
    a <- x %*% t(Linv)

    nd <- if (missing(newdata)) data.frame(1) else newdata

    ret <- lapply(1:J, function(j) {
        tmp <- object$models$models[[j]]
        prb <- tmp$todistr$p(a[,j])
        cf <- coef(tmp)
        cf[] <- object$models$parm(coef(object))[[j]]
        coef(tmp) <- cf
        predict(tmp, newdata = nd, type = "quantile", prob = prb)
    })
    
    ret <- do.call("cbind", ret)
    return(ret)
}

HDR <- function(object, level = .95, ...)
    UseMethod("HDR")

HDR.Mmlt <- function(object, level = .95, newdata, nsim = 1000L, K = 25, ...) {

    if (!missing(newdata)) {
        stopifnot(nrow(newdata) == 1)
    } else {
        newdata <- data.frame(1)
    }

    ### https://doi.org/10.2307/2684423 Section 3.2
    y <- simulate(object, newdata = newdata[rep(1, nsim),,drop = FALSE])
    y <- cbind(y, newdata)
    d <- predict(object, newdata = y, type = "density")

    ret <- do.call("expand.grid", lapply(object$models$models, 
                                         function(x) mkgrid(x, n = K)[[1L]]))
    colnames(ret) <- variable.names(object, response_only = TRUE)
    ret <- cbind(ret, newdata)
    ret$density <- predict(object, newdata = ret, type = "density")
    attr(ret, "cuts") <- quantile(d, prob = 1 - level)
    ret
}
