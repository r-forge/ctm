
mkgrid.cotram <- function(object, n, ...)
  mkgrid(object$count_response, n = n, ...)

as.mlt.cotram <- function(object) {
  cls <- which(class(object) == "mlt_fit")
  class(object) <- class(object)[-(1:(cls - 1))]
  object
}

logLik.cotram <- function(object, parm = coef(as.mlt(object), fixed = FALSE), newdata, ...){
  response <- variable.names(object, "response")
  
  if (!missing(newdata)) {
    ## check whether response is non-negative integer
    if (any(newdata[, response] < 0))
      stop("response is non-positive")
    if (!all(newdata[, response] %% 1 == 0))
      stop("response is non-integer")
    
    newdata[, response] <- as.integer(newdata[, response]) + as.integer(object$log_first)
    return(logLik(as.mlt(object), parm = parm, newdata = newdata, ...))
  }
  logLik(as.mlt(object), parm = parm, ...)
}
