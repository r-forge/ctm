
cotram <- function(formula, data, method = c("logit", "cloglog", "loglog", "probit"),
                   log_first = TRUE, prob = 0.9, subset, weights, offset, cluster,
                   na.action = na.omit, ...)
{
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action", "weights", "offset", "cluster"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(tram_data)
    td <- eval(mf, parent.frame())
    
    method <- match.arg(method)
    distribution <- c("logit" = "Logistic", "probit" = "Normal", 
                      "loglog" = "MaxExtrVal", "cloglog" = "MinExtrVal")
    distribution <- distribution[method]
    name <- c("logit" = "Odds", "loglog" = "Reverse Time Hazards",
              "cloglog" = "Hazards Cox")
    
    stopifnot(inherits(td$response, "response") || is.numeric(td$response))
    
    ## check that response is positive integer
    .check_count_var(td$response)
    # if (any(td$response < 0))
    #     stop("response is not a positive number")
    # if (!all(td$response %% 1 == 0))
    #     stop("response is not an integer number")
    
    y <- as.integer(td$response)
    
    ## y + 1 for log_first
    stopifnot(is.logical(log_first))
    plus_one <- as.integer(log_first)
    
    ## interval-censored count response for correct likelihood
    td$response <- .count_var(td$response, plus_one = plus_one)
    td$mf[, td$rname] <- .count_var(td$mf[, td$rname], plus_one = plus_one)
    
    # support & bounds
    support <- c(0, quantile(y, probs = prob)) + plus_one
    bounds <- c(-0.01, Inf) + plus_one
    
    ret <- tram(td, transformation = "smooth", distribution = distribution, 
                log_first = log_first, support = support, bounds = bounds, 
                negative = TRUE, ...)
    if (!inherits(ret, "mlt")) return(ret)
    
    ret$call <- match.call(expand.dots = TRUE)
    ret$log_first <- log_first
    ret$support <- support
    ret$bounds <- bounds
    if (method != "probit") {
        ret$tram <- paste(ifelse(is.null(td$terms$s), "", "(Stratified)"), 
                          "Discrete", name[method], "Count Transformation Model")
    } else {
        ret$tram <- paste(ifelse(is.null(td$terms$s), "", "(Stratified)"),
                          "Transformed Counts Probit Transformation Model")
    }
    
    ## <FIXME> return Surv object or count response? <\FIXME>
    ret$data[, td$rname] <- y
    ret$count_response <- numeric_var(td$rname, support = min(y):max(y))
    class(ret) <- c("cotram", class(ret))
    ret
}

.count_var <- function(y, plus_one) {
  ## count response as interval-censored object
  y <- as.integer(y)
  yleft <- y - 1L
  yleft[yleft < 0] <- -Inf
  
  Surv(yleft + plus_one, y + plus_one, type = "interval2")
}

.check_count_var <- function(y) {
  ## check that response is positive integer
  if (any(y < 0))
    stop("response is non-positive")
  if (!all(y %% 1 == 0))
    stop("response is non-integer")
}
