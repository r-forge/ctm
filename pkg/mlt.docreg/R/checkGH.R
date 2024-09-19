
checkGH <- function(object, tol = 1e-4) {

    nll <- function(parm) -logLik(object, parm = parm)

    ### check gradient  and hessian
    suppressWarnings(gr <- numDeriv::grad(nll, coef(object))) 
    s <- Gradient(object)
    ret <- isTRUE(all.equal(gr, s, check.attributes = FALSE, tol = tol))

    suppressWarnings(H1 <- numDeriv::hessian(nll, coef(object)))
    H2 <- Hessian(object)
    ret <- ret & isTRUE(all.equal(H1, H2, check.attributes = FALSE, tol = tol))
    return(ret)
}
