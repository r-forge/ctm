
### compute predicted values of (multiple) basis objects
predict.basis <- function(object, newdata, coef, 
                          dim = !is.data.frame(newdata), ...) {

    if (isTRUE(dim))
        dim <- sapply(newdata, NROW) 
    else if (is.logical(dim)) 
        dim <- NULL
    
    X <- model.matrix(object = object, data = newdata, dim = dim, ...)
    if (is.matrix(coef)) {
        if (is.null(dim) || nrow(coef) == nrow(X)) {
            lp <- rowSums(X * coef)
        } else {
            lp <- X %*% t(coef)
            stopifnot(max(abs(dim - dim(lp))) == 0)
            return(lp)
        }
    } else {
        lp <- c(X %*% coef)
    }
    if (is.null(dim)) return(lp)
    nd <- names(dim)
    ### <FIXME> essentially handle the length(dim) == 2 case
    if (any(nd %in% variable.names(object))) {
        nd <- nd[nd %in% variable.names(object)]
    } else {
        stopifnot(length(dim) == 2)
        nd <- names(dim)[2]
    }
    ### <FIXME>
    return(.const_array(dim, nd, lp))
}

nparm <- function(object)
    UseMethod("nparm")

nparm.basis <- function(object) {
    tmpdata <- mkgrid(as.vars(object), n = 1)
    tmpdata <- lapply(tmpdata, function(x) x[1])
    ncol(model.matrix(object, 
                      data = as.data.frame(tmpdata, 
                                           check.names = FALSE))) 
}

nparm.intercept_basis <- function(object)
    return(1)

nparm.box_bases <- function(object)
    prod(sapply(object, nparm))

nparm.cbind_bases <- function(object)
    sapply(object, nparm)

variable.names.basis <- function(object, ...)
    variable.names(as.vars(object))

variable.names.bases <- function(object, ...)
    unique(unlist(lapply(object, variable.names)))

mkgrid.basis <- function(object, n, ...)
    mkgrid(as.vars(object), n = n, ...)

mkgrid.bases <- mkgrid.basis

mkgrid.intercept_basis <- function(object, n, ...)
    return(c("(Intercept)" = 1L))

as.vars.basis <- function(object)
    attr(object, "variables")

as.vars.bases <- function(object) {
    vn <- variable.names(object)
    ret <- vector(mode = "list", length = length(vn))
    names(ret) <- vn
    .vars <- function(o) {
        if (inherits(o, "basis")) {
            v <- as.vars(o)
            if (inherits(v, "var")) v <- c(v)
            ret[variable.names(v)] <<- v
        } else 
            sapply(o, .vars)
    }
    ret2 <- .vars(object)
    do.call("c", ret)
}

bounds.basis <- function(object)
    bounds(as.vars(object))

bounds.bases <- function(object)
    bounds(as.vars(object))

as.basis <- function(object, ...)
    UseMethod("as.basis")

intercept <- function(object, ...)
    UseMethod("intercept")

intercept.default <- function(object, ...)
    attr(object, "intercept")

intercept.box_bases <- function(object, ...)
    any(sapply(object, intercept))

intercept.cbind_bases <- function(object, ...)
    any(sapply(object, intercept))
