
### marginally interpretable linear transformation models for clustered
### observations

mtram <- function(object, formula, data, 
                  grd = SparseGrid::createSparseGrid(type = "KPU", dimension = length(rt$cnms[[1]]), 
                                                     k = 10),
                  Hessian = FALSE, tol = .Machine$double.eps,
                  # standardise = FALSE,
                  ...) {
    
    standardise <- TRUE
    stopifnot(inherits(object, "mlt_fit"))
    
    bar.f <- lme4::findbars(formula)
    mf <- model.frame(lme4::subbars(formula), data = data)
    rt <- lme4::mkReTrms(bar.f, mf)
    
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
        PF <- function(z) qnorm(pmin(1 - tol, pmax(tol, P(z))))
        f <- object$todistr$d
        fprime <- object$todistr$dd
    }
    
    gr <- NULL
    
    ### catch constraint violations here
    .log <- function(x) 
        log(pmax(.Machine$double.eps, x))
    
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
            logdet <- 2 * determinant(L, logarithm = TRUE)$modulus
            
            # Sigma <- tcrossprod(LM)
            #            SigmaInv <- crossprod(Linv)
            D <- Dinv <- 1L
            if(standardise) {
                # D <- sqrt(diag(Sigma))
                D <- sqrt(rowSums(LM^2))
                Dinv <- 1/D
            }
            
            z <-  c(eY$Y %*% theta + offset)
            Dinvz <- Dinv * z
            PFz <- c(PF(Dinvz))
            DPFz <- c(D * PFz)
            
            if (NORMAL) {
                ret <- -0.5 * (logdet + sum((Linv %*% DPFz)^2) - sum(DPFz^2)) -
                    object$loglik(theta, weights = w)
            } else {
                ret <- -0.5 * (logdet + sum((Linv %*% DPFz)^2) - sum(PFz^2)) +
                    sum(.log(f(Dinvz) * c(eY$Yprime %*% theta)))
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
                        sum(c(fprime(Dinvz))/c(f(Dinvz)) * c(dDinv * z))
                }
                dtheta <- - crossprod(DPFz, SiID2) %*%
                    ((c(f(Dinvz))/dnorm(PFz)) * eY$Y) +
                    colSums((c(fprime(Dinvz))/c(f(Dinvz))) * (Dinv * eY$Y)) +
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
        
        ll <- function(parm) {
            theta <- parm[1:ncol(iY$Yleft)]
            gamma <- parm[-(1:ncol(iY$Yleft))]
            Lambdat@x[] <- mapping(gamma)
            lplower <- c(iY$Yleft %*% theta + offset)
            lplower[is.na(lplower)] <- -Inf
            lpupper <- c(iY$Yright %*% theta + offset)
            lpupper[is.na(lpupper)] <- Inf
            
            ## don't spend time on Matrix dispatch
            mLt <- t(as(Lambdat[wh, wh], "matrix"))
            ONE <- matrix(1, nrow = NCOL(mLt))
            
            ret <- sapply(1:length(idx), function(i) {
                V <- zt[[i]] %*% mLt  ### = U_i %*% Lambda(\varparm)
                i <- idx[[i]]
                if (standardise) {
                    sd <- c(sqrt((V^2) %*% ONE + 1)) ### D(\varparm)
                    zlower <- PF(lplower[i] / sd) * sd
                    zupper <- PF(lpupper[i] / sd) * sd
                } else {
                    zlower <- PF(lplower[i])
                    zupper <- PF(lpupper[i])
                    sd <- 1
                }
                .Marsaglia_1963(zlower, zupper, mean = 0, V = V, 
                                do_qnorm = FALSE, grd = grd)
            })
            return(-sum(.log(ret)))
        }
        X <- iY$Yleft
    }            
    
    ui <- attr(X, "constraint")$ui[, wf, drop = FALSE]
    ci <- attr(X, "constraint")$ci
    ui <- as(bdiag(ui, Diagonal(length(theta))), "matrix")
    ci <- c(ci, rt$lower)
    
    start <- c(coef(as.mlt(object), fixed = FALSE), theta)
    
    if (is.null(gr)) {
        opt <- alabama::auglag(par = start, fn = ll, 
                               hin = function(par) ui %*% par - ci, 
                               hin.jac = function(par) ui,
                               control.outer = list(trace = FALSE))[c("par", "value", "gradient")]
    } else {
        opt <- alabama::auglag(par = start, fn = ll, gr = gr,
                               hin = function(par) ui %*% par - ci, 
                               hin.jac = function(par) ui,
                               control.outer = list(trace = FALSE))[c("par", "value", "gradient")]
    }

    gamma <- opt$par[-(1:ncol(X))]
    names(opt$par)[-(1:ncol(X))] <- paste0("gamma", 1:length(gamma))
    Lambdat@x[] <- mapping(gamma)
    opt$G <- crossprod(Lambdat)[1:length(rt$cnms[[1]]),1:length(rt$cnms[[1]])]
    if (Hessian) opt$Hessian <- numDeriv::hessian(ll, opt$par)
    opt$loglik <- ll
    if(!is.null(gr)) opt$gr <- gr
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

.Marsaglia_1963 <- function(lower = rep(-Inf, nrow(sigma)), 
                            upper = rep(Inf, nrow(sigma)), 
                            mean = rep(0, nrow(sigma)), 
                            V = diag(2), 
                            grd = NULL,
                            do_qnorm = TRUE,
                            ...) {
    
    k <- nrow(V)
    l <- ncol(V)
    
    if (is.null(grd)) {
        stopifnot(do_qnorm)
        grd <- SparseGrid::createSparseGrid(type = "KPU", dimension = ncol(V), k = 10)
        ### Note: We expect t(nodes) below
        grd$nodes <- t(grd$nodes)
    }
    
    ### Note: The appendix describes the standardised version
    ### in order to be compliant with the notation in Genz & Bretz (2009).
    ### Standardisation of lower/upper and V AND division by diagonal
    ### elements in the Marsaglia formula cancel out and are thus
    ### left-out in the code
    
    lower <- lower - mean
    upper <- upper - mean
    
    if (k == 1) {
        VVt <- base::tcrossprod(V)
        sd <- sqrt(base::diag(VVt) + 1)
        return(pnorm(upper / sd) - pnorm(lower / sd))
    }
    
    ### y = qnorm(x)
    inner <- function(y) {
        Vy <- V %*% y
        ### this needs ~ 75% of the total runtime
        #ret <- pnorm(upper - Vy) - pnorm(lower - Vy)
        ### ~ 3x speed-up
        ### ret <- .Call("pnormMRS", c(upper - Vy)) - .Call("pnormMRS", c(lower - Vy))
        #if (nrow(Vy) == 1) return(ret)
        #ret <- matrix(pmax(.Machine$double.eps, ret), nrow = nrow(Vy),
        #              ncol = ncol(Vy))
        #exp(colSums(log(ret)))
        .Call("R_inner", upper - Vy, lower - Vy)
    }
    
    if (do_qnorm) grd$nodes <- qnorm(grd$nodes)
    ev <- inner(grd$nodes)
    c(value = sum(grd$weights * ev))
}
