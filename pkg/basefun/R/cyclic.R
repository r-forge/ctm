
cyclic_basis <- function(var, order = 3, frequency) {

    stopifnot(inherits(var, "numeric_var"))
    varname <- variable.names(var)
    support <- support(var)[[varname]]
    bounds <- bounds(var)[[varname]]

    basis <- function(data, deriv = 0L) {

        stopifnot(check(var, data))
        if (is.atomic(data)) {
            x <- data
        } else {
            if (is.null(varname)) varname <- colnames(data)[1]
            x <- data[[varname]]
        }
        if (frequency < bounds[1L] || frequency > bounds[2L])
            warning("bounded variable with unmatching frequency")

        X <- matrix(2 * pi * (1:order) / frequency, 
                    nrow = length(x), ncol = order, byrow = TRUE) * c(x)
        if (deriv == 0L)
            X <- cbind(sin(X), cos(X))
        if (deriv == 1L)
            X <- cbind(cos(X), -sin(X))
        if (deriv == 2L)
            X <- cbind(-sin(X), -cos(X))
        if (deriv == 3L)
            X <- cbind(-cos(X), sin(X))
        colnames(X) <- c(paste0("sin", 1:order), paste0("cos", 1:order))
        attr(X, "constraint") <- NULL
        attr(X, "Assign") <- matrix(varname, ncol = ncol(X))
        X
    }

    attr(basis, "variables") <- var
    attr(basis, "intercept") <- TRUE

    class(basis) <- c("cyclic_basis", "basis", class(basis))
    return(basis)
}

model.matrix.cyclic_basis <- function(object, data, deriv = 0L, ...)
    object(data, deriv = .deriv(variable.names(object), deriv))
