
mkgrid.cotram <- function(object, n, ...)
  mkgrid(object$count_response, n = n, ...)

# mkgrid.cotram <- function(object, n, smooth = FALSE, ...) {
#   ret <- mkgrid(object$count_response, ...)
#   y <- variable.names(object, "response")
#   if (smooth) ret[[y]] <-
#       seq(from = min(ret[[y]]), to = max(ret[[y]]), length.out = n)
#   ret
# }

as.mlt.cotram <- function(object) {
  class(object) <- class(object)[-which(class(object) == "cotram")]
  # if (object$log_first) warning("The model was fitted to 'response + 1L'.")
  as.mlt(object)
}

logLik.cotram <- function(object, parm = coef(as.mlt(object), fixed = FALSE), newdata, ...){
  response <- variable.names(object, "response")
  if (!missing(newdata)) {
    ## check that response is non-negative integer
    if (any(newdata[, response] < 0)) stop("response is non-positive")
    if (!all(newdata[, response] %% 1 == 0)) stop("response is non-integer")
    
    newdata[, response] <- R_count(y = newdata[, response], plus_one = as.integer(object$log_first))
    return(logLik(as.mlt(object), parm = parm, newdata = newdata, ...))
  }
  logLik(as.mlt(object), parm = parm, ...)
}

R_count <- function(y, plus_one) {
  ## count response as interval-censored object
  y <- as.integer(y)
  yleft <- y - 1L
  yleft[yleft < 0] <- -Inf
  
  Surv(yleft + plus_one, y + plus_one, type = "interval2")
}
