
### convert lp to PI and back, use logistic as default
PI.default <- function(object, prob, link = "logistic", ...) {

    FUN <- .lp2PI(link = link)

    if (missing(prob))
        return(FUN(object))

    object <- 1:999 / 50
    s <- spline(x = object, y = FUN(object), method = "hyman")
    wl5 <- (prob < .5 - .Machine$double.eps)
    wg5 <- (prob > .5 + .Machine$double.eps)
    ret <- numeric(length(prob))
    if (any(wl5))
        ret[wl5] <- - approx(x = s$y, y = s$x, xout = 1 - prob[wl5])$y
    if (any(wg5))
        ret[wg5] <- approx(x = s$y, y = s$x, xout = prob[wg5])$y
    ret
}
