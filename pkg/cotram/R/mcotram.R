
mcotram <- function(..., formula = ~ 1, data, conditional = FALSE, 
                    theta = NULL, fixed = NULL, scale = FALSE, 
                    optim = mmltoptim(), M = 1000,
                    dofit = TRUE, domargins = TRUE) {

    m <- list(...)
    stopifnot(all(sapply(m, function(x) inherits(x, "cotram"))))
    m <- lapply(m, as.mlt)

    J <- length(m)
    W <- t(ghalton(M, d = J - 1))

    m$formula <- formula
    m$data <- data
    m$conditional <- conditional
    m$theta <- theta
    m$fixed <- fixed
    m$scale <- scale
    m$optim <- optim
    m$args <- list(M = M, w = W)
    m$dofit <- dofit
    m$domargins <- domargins
    ret <- do.call("mmlt", m)
    class(ret) <- c("mcotram", class(ret))
    return(ret)
}

simulate.mcotram <- function(object, nsim = 1, seed = NULL, ...) {
    class(object) <- class(object)[-1L]
    ret <- as.double(simulate(object = object, nsim = nsim, seed = seed, ...))
    ret <- ceiling(ret)
    storage.mode(ret) <- "integer"
    plus_one <- sapply(object$models$models, function(x) x$log_first)
    if (any(plus_one)) {
        po <- matrix(plus_one, nrow = nrow(ret), ncol = ncol(ret), byrow = TRUE)
        ret <- pmax(ret - po, 0L)
    }
    return(ret)
}
