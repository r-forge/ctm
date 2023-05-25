
confband.cotram <- function(object, newdata,  level = 0.95, 
                            type = c("trafo", "distribution", "survivor", "cumhazard"), 
                            smooth = FALSE, q = NULL, K = 20, cheat = K, ...){
    stopifnot(!missing(newdata))
    stopifnot(is.data.frame(newdata))
    type <- match.arg(type)
    if (nrow(newdata) > 1) {
        ret <- lapply(1:nrow(newdata), function(i)
            confband(object = object, newdata = newdata[i,, drop = FALSE], q = q,
                     level = level, type = type, K = K, cheat = K, ...))
        return(ret)
    }
    stopifnot(nrow(newdata) == 1)
    
    y <- variable.names(object, "response")
    
    ## (q + 1) for log_first = TRUE
    plus_one <- as.integer(object$log_first)
    
    if (!is.null(q)) {
        if (any(q < 0)) stop("q is non-positive")
        if (!smooth && !all(q %% 1 == 0)) stop("q is non-integer")
        q <- q + plus_one
    } else {
        q <- mkgrid(object, n = K)[[y]] + plus_one
        
        if (smooth)
            q <- seq(from = min(q), to = max(q), length.out = K)
        
        if (!smooth && length(q) > K){
            cheat <- length(q)
            q <- floor(seq(from = min(q), to = max(q), length.out = K))
        }
    } 
    
    nd <- newdata[rep(1, length(q)),, drop = FALSE]
    nd[[y]] <- q
    
    X <- model.matrix(object$model$model, data = nd)
    
    ci <- confint(multcomp::glht(multcomp::parm(coef(as.mlt(object)), vcov(as.mlt(object))),
                                 linfct = X), ...)$confint
    
    ## use quantile obtained for K contrasts for larger number of contrasts
    if (cheat > K)
        return(confband(object = object, newdata = newdata, level = level,
                        type = type, K = cheat, calpha = attr(ci, "calpha")))
    
    if (type == "distribution") ci <- object$model$todistr$p(ci)
    if (type == "survivor") ci <- 1 - object$model$todistr$p(ci)
    if (type == "cumhazard") ci <- -log(1 - object$model$todistr$p(ci))
    
    ci <- cbind(q = q - plus_one, ci)
    return(ci)
}