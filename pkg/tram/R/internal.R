
asvar <- function(object, name, ...)
    UseMethod("asvar")

asvar.factor <- function(object, name, ...)  
    variables::factor_var(name, levels = levels(object))

asvar.ordered <- function(object, name, ...)  
    variables::ordered_var(name, levels = levels(object))

asvar.numeric <- function(object, name, prob = c(.1, .9), support = NULL, ...) {
    if (is.integer(object))
        return(variables::numeric_var(name, support = sort(unique(object))))  
    if (is.null(support)) support = quantile(object, prob = prob)
    variables::numeric_var(name, support = support, ...)
}

asvar.Surv <- function(object, name, prob = c(.1, .9), support = NULL, bounds = c(0, Inf), ...) {
    if (is.null(support)) {
        support <- quantile(survfit(y ~ 1, data = data.frame(y = object)), prob = prob)$quantile
        if (is.na(support[2])) 
            support[2] <- max(object[, 1])
    }
    variables::numeric_var(name, support = support, bounds = bounds, ...)
}

mkbasis <- function(yvar, transformation = c("discrete", "linear", "logarithmic", "smooth"), 
                    order = 6) {

  transformation <- match.arg(transformation)

  if (inherits(yvar, "numeric_var")) {
      return(switch(transformation, 
          "discrete" = stop("Discrete transformation not defined for numeric variable"),
          "linear" =
              polynomial_basis(yvar, coef = c(TRUE, TRUE), 
                               ui = matrix(c(0, 1), nrow = 1), ci = 0),
          "logarithmic" = 
              log_basis(yvar, ui = "increasing"),
          "smooth" = 
              Bernstein_basis(yvar, ui = "increasing", order = order)
          )
      )
  }

  return(as.basis(yvar))
}
