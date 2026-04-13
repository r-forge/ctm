
cotram <- function(formula, data, method = c("logit", "cloglog", "loglog", "probit"),
                   log_first = NULL, prob = NULL, order = NULL, subset, weights, offset, cluster,
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
    if (any(td$response < 0)) stop("response is non-positive")
    if (!all(td$response %% 1 == 0)) stop("response is non-integer")
    
    y <- as.integer(td$response)
    
    ### linear h only makes sense with log_first
    if (!is.null(order)) {
        if (order == 1L) log_first <- TRUE
    }

    if (isFALSE(log_first)) log_first <- NULL
    ## y + 1 for log_first
    plus_one <- isTRUE(log_first) + 0L
    
    ## interval-censored count response for correct likelihood
    td$response <- R(.R.count(td$response, plus_one = plus_one))
    td$mf[, td$rname] <- .R.count(td$mf[, td$rname], plus_one = plus_one)
    
    ## support
    support <- NULL
    if (!is.null(prob)) 
        support <- c(0, round(quantile(y, probs = prob)))[0:1 + length(prob)] + plus_one
    
    ## bounds
    bounds <- NULL
    if (!is.null(prob))
        bounds <- c(-0.01, Inf) + plus_one
    
    if (is.null(order)) {
        ret <- tram(td, transformation = "smooth", distribution = distribution, 
                    log_first = log_first, support = support, bounds = bounds, 
                    negative = TRUE, ...)
    } else {
        ret <- tram(td, transformation = "smooth", distribution = distribution, 
                    log_first = log_first, support = support, bounds = bounds, 
                    order = order,
                    negative = TRUE, ...)
    }

    if (!inherits(ret, "mlt")) return(ret)
    
    ret$call <- match.call(expand.dots = TRUE)
    ret$log_first <- isTRUE(log_first)
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
    
    if (ret$log_first) message(paste0("Model was fitted to log(", td$rname, " + 1)."))
    ret
}
