
.pmax <- function(x, y) {
    y[y < x] <- x
    return(y)
}

.log <- function(x) {
    ret <- log(.pmax(.Machine$double.eps, x))
    dim(ret) <- dim(x)
    return(ret)
}

### marginally interpretable linear transformation models for clustered
### observations

mtram <- function(object, formula, data, 
                  grd = SparseGrid::createSparseGrid(type = "KPU", dimension = length(rt$cnms[[1]]), 
                                                     k = 10),
                  tol = .Machine$double.eps,
                  optim = mltoptim(hessian = TRUE),
                  ...) {

    call <- match.call()
    
    stopifnot(inherits(object, "mlt_fit"))
    
    bar.f <- reformulas::findbars(formula)
    mf <- model.frame(reformulas::subbars(formula), data = data)
    rt <- reformulas::mkReTrms(bar.f, mf)
    
    ZtW <- rt$Zt
    Lambdat <- rt$Lambdat
    Lind <- rt$Lind
    mapping <- function(theta)
        theta[Lind]
    theta <- rt$theta
    
    eY <- get("eY", environment(object$loglik))
    iY <- get("iY", environment(object$loglik))
    fixed <- get("fixed", environment(object$loglik))
    offset <- get("offset", environment(object$loglik))
    wf <- 1:length(coef(as.mlt(object)))
    
    ### exact continuous obs
    if (!is.null(eY)) {
        tmp <- attr(eY$Y, "constraint")
        wf <- !colnames(eY$Y) %in% names(fixed)
        eY$Y <- eY$Y[, wf,drop = FALSE]
        attr(eY$Y, "constraint") <- tmp
    }
    
    ### discrete or censored
    if (!is.null(iY)) {
        tmp <- attr(iY$Yleft, "constraint")
        wf <- !colnames(iY$Yleft) %in% names(fixed)
        iY$Yleft <- iY$Yleft[, wf,drop = FALSE]
        iY$Yright <- iY$Yright[, wf,drop = FALSE]
        attr(iY$Yleft, "constraint") <- tmp
        attr(iY$Yright, "constraint") <- tmp
    }
    if ((length(eY$which) > 0) && (length(iY$which) > 0))
        stop("cannot deal with mixed censoring")
    
    w <- object$weights
    
    NORMAL <- FALSE
    if (object$todistr$name == "normal") {
        NORMAL <- TRUE
        PF <- function(z) z
    } else {
        P <- object$todistr$p
        PF <- function(z) qnorm(P(z, log.p = TRUE), log.p = TRUE) ### qnorm(pmin(1 - tol, pmax(tol, P(z))))
        f <- object$todistr$d
        ### fprime() / f()
        fpf <- object$todistr$dd2d
        dPF <- function(z)
            f(z) / dnorm(PF(z))
    }
    
    gr <- NULL
    
    ## continuous case
    if (length(eY$which) > 0) {
        L <- Cholesky(crossprod(Lambdat %*% ZtW), LDL = FALSE, Imult = 1)
        
        ll <- function(parm) {
            theta <- parm[1:ncol(eY$Y)]
            gamma <- parm[-(1:ncol(eY$Y))]
            
            Lambdat@x[] <- mapping(gamma)
            L <- update(L, t(Lambdat %*% ZtW), mult = 1)
            LM <- as(L, "Matrix")
            Linv <- solve(LM)
            logdet <- determinant(L, logarithm = TRUE, sqrt = FALSE)$modulus
            
            # Sigma <- tcrossprod(LM)
            #            SigmaInv <- crossprod(Linv)
            D <- Dinv <- 1L
            if (!NORMAL) {
                # D <- sqrt(diag(Sigma))
                D <- sqrt(rowSums(LM^2))
                Dinv <- 1/D
            }
            
            z <-  c(eY$Y %*% theta + offset)
            Dinvz <- Dinv * z
            PFz <- c(PF(Dinvz))
            DPFz <- c(D * PFz)
            
            if (NORMAL) {
                ret <- -0.5 * (logdet + sum((Linv %*% DPFz)^2) - sum(DPFz^2)) +
                    object$loglik(theta, weights = w)
            } else {
                ret <- -0.5 * (logdet + sum((Linv %*% DPFz)^2) - sum(PFz^2)) +
                    sum(f(Dinvz, log = TRUE) + .log(c(eY$Yprime %*% theta)))
            } 
            return(-ret)
        }
        
        X <- eY$Y
        
        if(NORMAL) { ## probit link
            gr <- function(parm) {
                theta <- parm[1:ncol(eY$Y)]
                gamma <- parm[-(1:ncol(eY$Y))]
                devLambda <- devSigma <- vector(mode = "list", 
                                                length = length(gamma))
                dgamma <- numeric(length(gamma))
                z <- c(eY$Y %*% theta + offset)
                Lambdat@x[] <- mapping(gamma)
                L <- update(L, t(Lambdat %*% ZtW), mult = 1)
                LM <- as(L, "Matrix")
                Linv <- solve(LM)
                # Sigma <- tcrossprod(LM)
                SigmaInv <- crossprod(Linv)
                LambdaInd <- t(Lambdat)
                LambdaInd@x[] <- 1:length(gamma)
                for (i in 1:length(gamma)) {
                    ### Wang & Merkle (2018, JSS) compute derivative of G!
                    ### We need derivative of Lambda!
                    dLtL <- (LambdaInd == i) %*% Lambdat
                    devLambda[[i]] <- dLtL + t(dLtL)
                    devSigma[[i]] <- crossprod(ZtW, devLambda[[i]] %*% ZtW)
                    t1 <- SigmaInv %*% devSigma[[i]]
                    t2 <- t1 %*% SigmaInv
                    dgamma[i] <- -.5 * (sum(diag(t1)) - crossprod(z, t2) %*% z)
                }
                dtheta <- - z %*% SigmaInv %*% eY$Y + 
                    colSums(eY$Yprime / c(eY$Yprime %*% theta))
                return(-c(c(as(dtheta, "matrix")), dgamma))
            }
        } else { ## no probit link
            gr <- function(parm) {
                theta <- parm[1:ncol(eY$Y)]
                gamma <- parm[-(1:ncol(eY$Y))]
                devLambda <- devSigma <- vector(mode = "list",
                                                length = length(gamma))
                
                Lambdat@x[] <- mapping(gamma)
                L <- update(L, t(Lambdat %*% ZtW), mult = 1)
                LM <- as(L, "Matrix")
                Linv <- solve(LM)
                # Sigma <- tcrossprod(LM)
                SigmaInv <- crossprod(Linv)
                D2 <- rowSums(LM^2) # diag(Sigma)
                D <- sqrt(D2)
                Dinv <- 1/D
                Dinv2 <- 1/D2
                LambdaInd <- t(Lambdat)
                LambdaInd@x[] <- 1:length(gamma)
                
                dgamma <- numeric(length(gamma))
                z <-  c(eY$Y %*% theta + offset)
                Dinvz <- Dinv * z
                PFz <- c(PF(Dinvz))
                DPFz <- c(D * PFz)
                
                for (i in 1:length(gamma)) {
                    ### Wang & Merkle (2018, JSS) compute derivative of G!
                    ### We need derivative of Lambda!
                    dLtL <- (LambdaInd == i) %*% Lambdat
                    devLambda[[i]] <- dLtL + t(dLtL)
                    devSigma[[i]] <- crossprod(ZtW, devLambda[[i]] %*% ZtW)
                    t1 <- SigmaInv %*% devSigma[[i]]
                    t2 <- t1 %*% SigmaInv
                    dDPFz <- .5 * Dinv2 * diag(devSigma[[i]]) *
                        (DPFz - (c(f(Dinvz)) / dnorm(PFz)) * z)
                    dDinv <- -.5 * D^(-3) * diag(devSigma[[i]])
                    dDinv2 <- - D2^(-2) * diag(devSigma[[i]])
                    t3 <- t2
                    diag(t3) <- diag(t3) + dDinv2
                    SiID2 <- SigmaInv
                    diag(SiID2) <- diag(SiID2) - Dinv2
                    dgamma[i] <- -.5 * (sum(diag(t1)) +
                                            crossprod(dDPFz, SiID2) %*% DPFz -
                                            crossprod(DPFz, t3) %*% DPFz +
                                            crossprod(DPFz, SiID2) %*% dDPFz) +
                        sum(c(fpf(Dinvz)) * c(dDinv * z))
                }
                dtheta <- - crossprod(DPFz, SiID2) %*%
                    ((c(f(Dinvz))/dnorm(PFz)) * eY$Y) +
                    colSums(c(fpf(Dinvz)) * (Dinv * eY$Y)) +
                    colSums(eY$Yprime / c(eY$Yprime %*% theta))
                
                return(-c(c(as(dtheta, "matrix")), dgamma))
            }
        }
    } else { ## censored and discrete case
        if (length(rt$flist) != 1L)
            stop("only one grouping factor allowed for discrete case")
        grp <- rt$flist[[1]]
        idx <- split(1:length(grp), grp)
        wh <- 1:length(rt$cnms[[1]])
        ### .Marsaglia_1963 expects t(nodes) !!!
        grd$nodes <- t(qnorm(grd$nodes))
        ## don't spend time on Matrix dispatch
        mZtW <- as(ZtW, "matrix")
        zt <- lapply(idx, function(i) {
            z <- mZtW[,i,drop = FALSE]
            t(z[base::rowSums(abs(z)) > 0,,drop = FALSE])
        })

        ### <FIXME> scale parameters if necessary </FIXME>
        ll <- function(parm) {
            theta <- parm[1:ncol(iY$Yleft)]
            gamma <- parm[-(1:ncol(iY$Yleft))]
            Lambdat@x[] <- mapping(gamma)
            lplower <- c(iY$Yleft %*% theta + offset)
            lplower[!is.finite(lplower)] <- -Inf
            lpupper <- c(iY$Yright %*% theta + offset)
            lpupper[!is.finite(lpupper)] <- Inf
            
            ## don't spend time on Matrix dispatch
            mLt <- t(as(Lambdat[wh, wh], "matrix"))
            
            ret <- sapply(1:length(idx), function(i) {
                V <- zt[[i]] %*% mLt  ### = U_i %*% Lambda(\varparm)
                i <- idx[[i]]
                if (!NORMAL) {
                    sd <- c(sqrt(rowSums(V^2) + 1)) ### D(\varparm)
                    zlower <- PF(lplower[i] / sd) * sd
                    zupper <- PF(lpupper[i] / sd) * sd
                } else {
                    zlower <- lplower[i]
                    zupper <- lpupper[i]
                }
                lpRR(lower = zlower, upper = zupper, mean = 0, B = V, 
                     Z = grd$nodes, weights = grd$weights, log.p = TRUE)
            })
            return(-sum(ret))
        }
        X <- iY$Yleft

        if (NORMAL) {
            gr <- function(parm) {
                theta <- parm[1:ncol(iY$Yleft)]
                gamma <- parm[-(1:ncol(iY$Yleft))]
                Lambdat@x[] <- mapping(gamma)
                lplower <- c(iY$Yleft %*% theta + offset)
                lplower[!is.finite(lplower)] <- -Inf
                lpupper <- c(iY$Yright %*% theta + offset)
                lpupper[!is.finite(lpupper)] <- Inf
            
                ## don't spend time on Matrix dispatch
                mLt <- t(as(Lambdat[wh, wh], "matrix"))
            
                ret <- lapply(1:length(idx), function(i) {
                    V <- (B <- zt[[i]]) %*% mLt  ### = U_i %*% Lambda(\varparm)
                    i <- idx[[i]]
                    zlower <- PF(lplower[i])
                    zupper <- PF(lpupper[i])
                    ret <- slpRR(lower = zlower, upper = zupper, mean = 0, B = V, 
                                 Z = grd$nodes, weights = grd$weights, log.p = TRUE)
                    dtheta <- colSums(ret$lower * iY$Yleft[i,,drop = FALSE], na.rm = TRUE) + 
                              colSums(ret$upper * iY$Yright[i,,drop = FALSE], na.rm = TRUE)
                    K <- ncol(V)
                    ind <- matrix(1:(K^2), nrow = K, byrow = TRUE)
                    dgamma <- as.vector(t(ret$B) %*% B)[ind[lower.tri(ind, diag = TRUE)]]
                    return(c(dtheta, dgamma))
                })
                return(-Reduce("+", ret))
            }
        } else {
            dsd <- function(gamma, B) {
                L <- diag(ncol(B))
                L[lower.tri(L, diag = TRUE)] <- gamma
                sd <- sqrt(rowSums((A <- B %*% L)^2) + 1)
                ret <- lapply(1:nrow(B), function(i)
                    as.vector(A[i,] %*% t(B[i,])))
                ret <- do.call("rbind", ret)
                ret / sd
            }

            gr <- function(parm) {
                theta <- parm[1:ncol(iY$Yleft)]
                gamma <- parm[-(1:ncol(iY$Yleft))]
                Lambdat@x[] <- mapping(gamma)
                lplower <- c(iY$Yleft %*% theta + offset)
                lplower[!is.finite(lplower)] <- -Inf
                lpupper <- c(iY$Yright %*% theta + offset)
                lpupper[!is.finite(lpupper)] <- Inf
            
                ## don't spend time on Matrix dispatch
                mLt <- t(as(Lambdat[wh, wh], "matrix"))
            
                ret <- lapply(1:length(idx), function(i) {
                    V <- (B <- zt[[i]]) %*% mLt  ### = U_i %*% Lambda(\varparm)
                    i <- idx[[i]]
                    sd <- c(sqrt(rowSums(V^2) + 1)) ### D(\varparm)
                    zlower <- PF(lsd <- lplower[i] / sd) * sd
                    zupper <- PF(usd <- lpupper[i] / sd) * sd
                    ret <- slpRR(lower = zlower, upper = zupper, mean = 0, B = V, 
                                 Z = grd$nodes, weights = grd$weights, log.p = TRUE)
                    dtheta <- colSums(ret$lower * dPF(lsd) * iY$Yleft[i,,drop = FALSE], 
                                      na.rm = TRUE) + 
                              colSums(ret$upper * dPF(usd) * iY$Yright[i,,drop = FALSE], 
                                      na.rm = TRUE)
                    K <- ncol(V)
                    dsdg <- dsd(gamma, B = B)
                    dgamma <- colSums((ret$lower * (dPF(lsd) * (-lsd / sd) * sd +
                                                    PF(lsd))) * dsdg, na.rm = TRUE) + 
                              colSums((ret$upper * (dPF(usd) * (-usd / sd) * sd +
                                                    PF(usd))) * dsdg, na.rm = TRUE)
                    dgamma <- dgamma + as.vector(t(ret$B) %*% B)
                    ind <- matrix(1:(K^2), nrow = K, byrow = TRUE)
                    idx <- ind[lower.tri(ind, diag = TRUE)]
                    return(c(dtheta, dgamma[idx]))
                })
                return(-Reduce("+", ret))
            }
        }
    }            
    
    ui <- attr(X, "constraint")$ui[, wf, drop = FALSE]
    ci <- attr(X, "constraint")$ci
    ui <- as(bdiag(ui, Diagonal(length(theta))), "matrix")
    ci <- c(ci, rt$lower)
    ui <- as(ui[is.finite(ci),,drop = FALSE], "matrix")
    ci <- ci[is.finite(ci)]
    
    start <- c(coef(as.mlt(object), fixed = FALSE), theta)
    
    for (i in 1:length(optim)) {
        opt <- optim[[i]](theta = start, 
                          f = ll,
                          g = gr,
                          ui = ui, ci = ci)
        if (opt$convergence == 0) break()
    }
    if (opt$convergence != 0)
        warning("Optimisation did not converge")

    gamma <- opt$par[-(1:ncol(X))]
    names(opt$par)[-(1:ncol(X))] <- paste0("gamma", 1:length(gamma))
    Lambdat@x[] <- mapping(gamma)
    opt$G <- as(crossprod(Lambdat)[1:length(rt$cnms[[1]]),
                                   1:length(rt$cnms[[1]])], "matrix")
    opt$loglik <- ll
    if(!is.null(gr)) opt$gr <- gr
    opt$call <- call
    class(opt) <- "mtram"
    opt
}


logLik.mtram <- function(object, parm = NULL, ...) {
    if (!is.null(parm)) {
        ret <- -c(object$loglik(parm))
    } else {
        ret <- -c(object$value)
    }
    attr(ret, "df") <- length(coef(object))
    class(ret) <- "logLik"
    ret
}

coef.mtram <- function(object, ...)
    object$par

Hessian.mtram <- function(object, ...) {
    H <- object$optim_hessian
    return(H)
}

### FIXME: one would want sd for _marginal_ parameters,
### that is, coef(object) / sqrt(1 + coef(object)["gamma"]^2) in the
### simplest case of repeated measurements
vcov.mtram <- function(object, ...) {
    class(object) <- c("mtram", "mlt")
    ret <- mlt:::vcov.mlt(object)
    colnames(ret) <- rownames(ret) <- names(coef(object))
    ret <- (ret + t(ret)) / 2
    return(ret)
}
