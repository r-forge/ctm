
.add_confband <- function(object, fill = "lightgrey", lwd = 1, col = "black") {

    q <- object[, "q"]
    lwr <- object[, "lwr"]
    upr <-  object[, "upr"]
    polygon(c(q, rev(q)), c(lwr, rev(upr)),
            border = NA, col = fill) 
    lines(object[, "q"], object[, "Estimate"], col = col, lwd = lwd)
}

plot.tram <- function(x, newdata = model.frame(x), 
    which = c("QQ-PIT", "baseline only", "distribution"), 
    confidence = c("none", "interval", "band"), level = 0.95, 
    K = 50, cheat = K, col = "black", fill = "lightgrey", lwd = 1, ...) {

    which <- match.arg(which)
    object <- as.mlt(x)
    y <- newdata[[variable.names(x, "response")]]

    censored <- inherits(y, "Surv") || inherits(y, "response")

    if (which != "distribution" && censored)
        stop("Cannot compute in-sample ", which, " for censored responses")
        
    cb <- NULL
    confidence <- match.arg(confidence)
    calpha <- switch(confidence, "none" = NULL,
                                 "interval" = univariate_calpha(),
                                 "band" = adjusted_calpha())
    confidence <- confidence != "none"

    ret <- switch(which, 
        "QQ-PIT" = {
            U <- predict(object, newdata = newdata, type = "distribution")
            qqplot(U, 1:length(U) / length(U), col = col, ...)
            qqline(U, distribution = qunif)
        },
        "baseline only" = {
            scf <- object$shiftcoef
            if (length(scf) > 0) {
                mobj <- as.mlt(object)
                cf <- coef(mobj)
                cf[scf] <- 0
                coef(mobj) <- cf
                plot(mobj, newdata = newdata, type = "trafo", col = col, lwd = lwd, ...)
                if (confidence)
                    cb <- confband(mobj, newdata = newdata, type = "trafo", 
                                   calpha = calpha, level = level, K = K, 
                                   cheat = cheat)
            } else {
                plot(object, newdata = newdata, type = "trafo", col = col, lwd = lwd, 
                     K = K, ...)
                if (confidence) 
                    cb <- confband(object, newdata = newdata, type = "trafo", 
                                   calpha = calpha, level = level, K = K, 
                                   cheat = cheat)
            }
        },
        "distribution" = {
            plot(object, newdata = newdata, col = col, lwd = lwd, K = K, ...)
            type <- list(...)$type
            if (is.null(type)) type <- "trafo"
            if (confidence)
                cb <- confband(object, newdata = newdata, calpha = calpha, 
                               level = level, K = K, cheat = cheat, type = type)
        })

    if (confidence) {
        if (length(fill) != NROW(newdata)) 
            fill <- rep(fill, length.out = NROW(newdata))
        if (length(col) != NROW(newdata)) 
            col <- rep(col, length.out = NROW(newdata))

        if (is.matrix(cb)) {
            .add_confband(cb, fill = fill[1], col = col[1], lwd = lwd) 
        } else {
            out <- lapply(1:length(cb), function(i) 
                .add_confband(cb[[i]], fill = fill[i], col = col[i], lwd = lwd))
        }
    }
}

plot.ROCtram <- function(x, lty = 1:ncol(x), 
    col = "black", fill = "lightgrey", lwd = 1, ...) {

    prob <- attr(x, "prob")
    plot(0, 1, xlim = c(0, 1), ylim = c(0, 1), type = "n", 
         xlab = "1 - Specificity", ylab = "Sensitivity", ...)
    abline(a = 0, b = 1, col = "lightgrey")
    if (length(fill) != ncol(x)) 
        fill <- rep(fill, length.out = ncol(x))
    if (length(col) != ncol(x)) 
        col <- rep(col, length.out = ncol(x))
    if (length(lty) != ncol(x)) 
        lty <- rep(lty, length.out = ncol(x))

    for (i in 1:ncol(x))
        lines(prob, x[, i], lty = lty[i])
    cb <- attr(x, "conf.band")
    if (!is.null(cb)) {
        for (i in 1:ncol(cb$lwr))
            .add_confband(cbind(q = prob, Estimate = x[,i], 
                                lwr = cb$lwr[, i], upr = cb$upr[, i]),
                          fill = fill[i], col = col[i], ...)
    }
}
