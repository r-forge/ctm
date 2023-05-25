
mcotram <- function(..., formula = ~ 1, data, conditional = FALSE, 
                    theta = NULL, optim = mmltoptim(), M = 1000,
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
    m$optim <- optim
    m$args <- list(M = M, w = W)
    m$dofit <- dofit
    m$domargins <- domargins
    return(do.call("mmlt", m))
}
