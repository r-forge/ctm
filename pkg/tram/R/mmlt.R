
Mmlt <- function(..., formula = ~ 1, data, conditional = FALSE, 
                 theta = NULL, fixed = NULL, scaleparm = FALSE,
                 optim = mltoptim(hessian = TRUE),
                 args = list(seed = 1, type = c("MC", "ghalton"), M = 1000), 
                 fit = c("jointML", "pseudo", "ACS", "sequential", "none"),
                 ACSiter = 2)
{

    call <- match.call()  
    fit <- match.arg(fit)

    if (fit %in% c("none", "jointML", "pseudo")) {
        ret <- switch(fit, 
            "none" = mmlt(..., formula = formula, data = data, conditional = conditional, 
                          theta = theta, fixed = fixed, scaleparm = scaleparm,
                          optim = optim, args = args,
                          dofit = FALSE),
            "jointML" = mmlt(..., formula = formula, data = data, conditional = conditional, 
                          theta = theta, fixed = fixed, scaleparm = scaleparm,
                          optim = optim, args = args),
            "pseudo" = mmlt(..., formula = formula, data = data, conditional = conditional, 
                          theta = theta, fixed = fixed, scaleparm = scaleparm,
                          optim = optim, args = args, domargins = FALSE))
        ret$call <- call
        class(ret) <- c("Mmlt", class(ret))
        return(ret)
    } 

    if (fit == "sequential") {
        m <- list(...)
        if (!is.null(theta))
            stop("Cannot perform sequential fit with starting values")
        if (!is.null(fixed))
            stop("Cannot perform sequential fit with fixed values")
        mj <- as.mlt(m[[1]])
        for (j in 2:length(m)) {
            mc <- m[1:j]
            mc$formula <- formula
            mc$data <- data
            mc$conditional <- conditional
            mc$scale <- scale
            mc$optim <- optim
            mc$args <- args
            mc$dofit <- FALSE
            ret <- do.call("mmlt", mc) ### only get names
            cfj <- coef(mj, fixed = TRUE)
            if (j == 2)
                names(cfj) <- ret$parnames[1:length(cfj)]
            ### fix coefs and estimate jth row of
            ### lambda and jth marginal parameters only
            mc$fixed <- cfj
            mc$dofit <- TRUE
            mj <- do.call("mmlt", mc)
        }
        mj$call <- call
        class(mj) <- c("Mmlt", class(mj))
        return(mj)
    }

    ### ACS: start with Lambda parameters, keep margins fix
    ret <- mmlt(..., formula = formula, data = data, conditional = conditional, 
               theta = theta, fixed = fixed, scaleparm = scaleparm,
               optim = optim, args = args, domargins = FALSE)

    cf <- coef(ret)
    P <- length(cf)
    pn <- names(cf)
    cf <- ret$parm(cf)
    nL <- length(cf[[length(cf)]])
    nM <- P - nL
    mn <- pn[1:nM]
    Ln <- pn[-(1:nM)]
    for (i in 1:ACSiter) {
        ### update marginal parameters given Lambda
        theta <- coef(ret)[mn]  ### marginal parameters
        fixed <- coef(ret)[Ln]  ### Lambda parameters
        ret <- mmlt(..., formula = formula, data = data, conditional = conditional, 
                    theta = theta, fixed = fixed, scaleparm = scaleparm,
                    optim = optim, args = args)
        ### update Lambda given marginal parameters
        theta <- coef(ret)[Ln]	### Lambda parameters
        fixed <- coef(ret)[mn]  ### marginal parameters
        ret <- mmlt(..., formula = formula, data = data, conditional = conditional, 
                    theta = theta, fixed = fixed, scaleparm = scaleparm,
                    optim = optim, args = args)
    }
    ret$call <- call
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
    y <- simulate(object, newdata = newdata[rep_len(1, nsim),,drop = FALSE])
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
