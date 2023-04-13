
### catch constraint violations here
.log <- function(x) {
  ret <- log(pmax(.Machine$double.eps, x))
  dim(ret) <- dim(x)
  ret
}

.chol <- function(Lambda) {
    chol <- solve(Lambda)
    CCt <- Tcrossprod(chol, diag_only = TRUE)
    ret <- Dchol(chol, D = 1 / sqrt(CCt))
    return(ret)
}

.magic <- function(Lambda, Dchol, N) {

    J <- dim(Lambda)[2L]

    ### vectrick needs byrow = FALSE, so do it here once
    Lambda <- ltMatrices(Lambda, byrow = FALSE)

    chol <- solve(Lambda)
    CCt <- Tcrossprod(chol, diag_only = TRUE)
    DC <- Dchol(chol, D = Dinv <- 1 / sqrt(CCt))
    SDC <- solve(DC)

    IDX <- t(M <- matrix(1:J^2, nrow = J, ncol = J))
    i <- cumsum(c(1, rep(J + 1, J - 1)))
    ID <- diagonals(as.integer(J), byrow = attr(Lambda, "byrow"))
    if (dim(ID)[1L] != dim(chol)[1L])
        ID <- ID[rep(1, dim(chol)[1L]),]

    if (inherits(Dchol, "ltMatrices")) {
        T1 <- matrix(as.array(Dchol), nrow = dim(Dchol)[2L]^2)
    } else {
        T1 <- Dchol
    }

    B <- vectrick(ID, T1, chol)
    B[i,] <- B[i,] * (-.5) * c(CCt)^(-3/2)
    B[-i,] <- 0

    Dtmp <- Dchol(ID, D = Dinv)

    ret <- vectrick(ID, B, chol, transpose = c(TRUE, FALSE)) +
           vectrick(chol, B, ID)[IDX,] +
           vectrick(Dtmp, T1, ID)

    ### this means: ret <- - vectrick(chol, ret, chol)
    ret <- - vectrick(chol, ret)
    ret <- ltMatrices(ret[M[lower.tri(M)],,drop = FALSE],
                      byrow = FALSE, diag = FALSE)
    ret <- ltMatrices(ret, diag = FALSE, byrow = TRUE)
    diagonals(ret) <- 0
    ret
}

.ll <- function(dim, scale = TRUE, args = list()) {

    if (length(dim) == 1L)
        dim <- c(dim, 0L)

    cJ <- dim[1L]
    dJ <- dim[2L]

    cll <- function(obs, Lambda) {

        stopifnot(!attr(Lambda, "diag"))

        if (!scale)
            return(dmvnorm(x = t(obs), invchol = Lambda, log = TRUE))

        chol <- .chol(Lambda)
        return(dmvnorm(x = t(obs), chol = chol, log = TRUE))
    }

    csc <- function(obs, Lambda, magic = TRUE) {

        # stopifnot(!attr(Lambda, "diag"))
        stopifnot(attr(Lambda, "byrow"))

        N <- ncol(obs)

        if (!scale) {
            ret <- sldmvnorm(x = t(obs), invchol = Lambda)
            names(ret)[names(ret) == "invchol"] <- "Lambda"
            names(ret)[names(ret) == "x"] <- "obs"
            return(ret)
        }

        if (!magic) {
            ret <- sldmvnorm(x = t(obs), invchol = invcholD(Lambda))
            names(ret)[names(ret) == "invchol"] <- "Lambda"
            names(ret)[names(ret) == "x"] <- "obs"
            return(ret)
        }

        chol <- .chol(Lambda)
        ret <- sldmvnorm(x = t(obs), chol = chol)
        dobs <- ret$x
        ret <- .magic(Lambda, ret$chol, N)
        return(list(Lambda = ret, obs = dobs))
    }

    if (!dJ) return(list(logLik = cll, score = csc))

    dll <- function(lower, upper, Lambda, center = NULL) {

        # stopifnot(!attr(Lambda, "diag"))

        a <- args
        a$center <- center
        a$mean <- 0
        a$lower <- lower
        a$upper <- upper
        a$logLik <- FALSE
        if (!scale || cJ > 0) {
            a$invchol <- Lambda
        } else {
            a$chol <- Dchol(solve(Lambda))
        }
        do.call("lpmvnorm", a)
    }

    dsc <- function(lower, upper, Lambda, center = NULL) {

        # stopifnot(!attr(Lambda, "diag"))
        stopifnot(attr(Lambda, "byrow"))

        a <- args
        a$center <- center
        a$mean <- 0
        a$lower <- lower
        a$upper <- upper
        if (!scale || cJ > 0) {
            a$invchol <- Lambda
            a$logLik <- TRUE
            ret <- do.call("slpmvnorm", a)
            names(ret)[names(ret) == "invchol"] <- "Lambda"
            return(ret)
        }

        a$chol <- .chol(Lambda)        
        ret <- do.call("slpmvnorm", a)
        smean <- ret$mean
        slower <- ret$lower
        supper <- ret$upper

        ret <- .magic(Lambda, ret$chol, ncol(lower))

        ret <- list(Lambda = ret, mean = smean, lower = slower, upper = supper)
        return(ret)
    }

    if (!cJ) return(list(logLik = dll, score = dsc))
        
    ll <- function(obs, lower, upper, Lambda) {

        md <- marg_mvnorm(invchol = Lambda, which = 1:cJ)
        ret <- cll(obs = obs, Lambda = md$invchol)

        if (scale) Lambda <- invcholD(Lambda)

        cd <- cond_mvnorm(invchol = Lambda, which_given = 1:cJ, given = obs, center = TRUE)
        ret <- ret + dll(lower = lower, upper = upper, Lambda = cd$invchol, 
                         center = cd$center)
        return(ret)
    }

    sc <- function(obs, lower, upper, Lambda) {

        if (scale) {
            L1 <- Lambda
            Lambda <- invcholD(Lambda)
        }
        md <- marg_mvnorm(invchol = Lambda, which = 1:cJ)
        cs <- csc(obs = obs, Lambda = md$invchol, magic = FALSE)

        cd <- cond_mvnorm(invchol = Lambda, which_given = 1:cJ, given = obs, center = TRUE)
        ds <- dsc(lower = lower, upper = upper, center = cd$center, Lambda = cd$invchol)

        tmp0 <- solve(cd$invchol, ds$mean, transpose = TRUE)
        tmp <- -tmp0[rep(1:dJ, each = cJ),,drop = FALSE] * obs[rep(1:cJ, dJ),,drop = FALSE]
        # tmp1 <- do.call("cbind", lapply(1:ncol(obs), function(i) -kronecker(tmp0[,i], obs[,i])))
        # stopifnot(isTRUE(all.equal(tmp, tmp1, check.attributes = FALSE)))

        J <- cJ + dJ
        Jp <- J * (J + c(-1, 1)[scale + 1L]) / 2
        M <- as.array(ltMatrices(1:Jp, diag = scale, byrow = TRUE))[,,1]
        ret <- matrix(0, nrow = Jp, ncol = ncol(obs))
        M1 <- M[1:cJ, 1:cJ]
        idx <- t(M1)[upper.tri(M1, diag = scale)]
        ret[idx,] <- Lower_tri(cs$Lambda, diag = scale)

        idx <- c(t(M[-(1:cJ), 1:cJ]))
        ret[idx,] <- tmp

        M3 <- M[-(1:cJ), -(1:cJ)]
        idx <- t(M3)[upper.tri(M3, diag = scale)]
        ret[idx,] <- Lower_tri(ds$Lambda, diag = scale)

        ret <- ltMatrices(ret, diag = scale, byrow = TRUE)

        if (scale) {
            k1 <- matrix(as.array(ret), nrow = J^2)
            ## this means: T1 <- -vectrick(Lambda, k1, Lambda)
            T1 <- -vectrick(Lambda, k1)
            ret <- .magic(L1, T1, ncol(obs))
        } else {
            diagonals(ret) <- 0
        }

        ### post differentiate mean 
        aL <- as.array(Lambda)[-(1:cJ), 1:cJ,]
        lst <- tmp0[rep(1:dJ, cJ),,drop = FALSE]
        dobs <- -margin.table(aL * array(lst, dim = dim(aL)), 2:3)

        ret <- c(list(Lambda = ret, obs = cs$obs + dobs), 
                 ds[c("lower", "upper")])
        return(ret)
    }

    return(list(logLik = ll, score = sc))
}

.models <- function(...) {

    m <- lapply(list(...), function(x) as.mlt(x))
    nm <- abbreviate(sapply(m, function(x) x$model$response), 4)
    J <- length(m)
    Jp <- J * (J - 1) / 2
    normal <- sapply(m, function(x) x$todistr$name == "normal")
  
    w <- lapply(m, weights)
    out <- lapply(w, function(x) stopifnot(isTRUE(all.equal(x, w[[1]]))))
    w <- w[[1L]]
    if (isTRUE(all.equal(unique(w), 1))) w <- FALSE
  
    mm <- lapply(m, function(mod) {
      eY <- get("eY", environment(mod$parm))
      iY <- get("iY", environment(mod$parm))
      list(eY = eY, iY = iY)
    })

    cmod <- sapply(mm, function(x) !is.null(x$eY))  
    dmod <- sapply(mm, function(x) !is.null(x$iY))  
    stopifnot(all(xor(cmod, dmod)))
    ### continuous models first
    stopifnot(all(diff(cmod) <= 0))
    stopifnot(all(diff(dmod) >= 0))

    nobs <- lapply(m, nobs)
    stopifnot(isTRUE(do.call("all.equal", nobs)))
    nobs <- nobs[[1L]]

    P <- sapply(m, function(x) length(coef(x)))
    fpar <- factor(rep(1:J, P))

    parm <- function(par) {
        mpar <- par[1:sum(P)]
        split(mpar, fpar)
    }

    constr <- lapply(mm, function(m) {
        if (is.null(m$eY)) return(attr(m$iY$Yleft, "constraint"))
        return(attr(m$eY$Y, "constraint"))
    })

    ui <- do.call("bdiag", lapply(constr, function(x) x$ui))
    ci <- do.call("c", lapply(constr, function(x) x$ci))
    ui <- as(ui[is.finite(ci),,drop = FALSE], "matrix")
    ci <- ci[is.finite(ci)]

    mf <- lapply(1:J, function(j) {
        mf <- m[[j]]$data ###model.frame(m[[j]])
        if (cmod[j]) return(mf)
        yl <- m[[j]]$response$cleft
        yr <- m[[j]]$response$cright
        rp <- m[[j]]$model$response
        ml <- mr <- mf
        ml[[rp]] <- yl
        mr[[rp]] <- yr
        return(list(left = ml, right = mr))
    })

    ### needs newdata for predict
    nn <- sapply(1:J, function(j) {
        !is.null(m[[j]]$fixed) ||
        !isTRUE(all.equal(unique(m[[j]]$offset), 0)) ||
        m[[j]]$model$scale_shift
    })

    type <- lapply(1:J, function(j)
        mlt:::.type_of_response(m[[j]]$response))

    return(list(models = m, mf = mf, cont = cmod, type = type, normal = normal, 
                nobs = nobs, weights = w, nparm = P, parm = parm, 
                ui = ui, ci = ci, mm = mm, names = nm, nn = nn))
}

.model_matrix <- function(models, j = 1, newdata = NULL, prime = FALSE) {

    if (is.null(newdata)) {
        if (models$cont[j]) {
            if (prime) return(models$mm[[j]]$eY$Yprime)
            return(models$mm[[j]]$eY$Y)
        }
        stopifnot(!prime)
        Yleft <- models$mm[[j]]$iY$Yleft
        Yright <- models$mm[[j]]$iY$Yright
        Yleft[!is.finite(Yleft[,1]),] <- 0
        Yright[!is.finite(Yright[,1]),] <- 0
        return(list(Yleft = Yleft, Yright = Yright))
    }

    if (models$cont[j]) {
        if (prime) {
            drv <- 1L
            names(drv) <- models$models[[j]]$model$response
            return(model.matrix(models$models[[j]]$model, data = newdata, deriv = drv))
        } else {
            return(model.matrix(models$models[[j]]$model, data = newdata))
        }
    }
    stop("not yet implemented")
}

.mget <- function(models, j = 1, parm, newdata = NULL,
                  what = c("trafo", "dtrafo", "z", "zleft", 
                           "dzleft", "zright", "dzright", "zprime", 
                           "mm", "mmprime", "estfun"), ...) {

    what <- match.arg(what)

    if (length(j) > 1) {
        ret <- lapply(j, .mget, models = models, parm = parm, 
                      newdata = newdata, what = what)
        return(ret)
    }

    prm <- models$parm(parm)[[j]]
    tmp <- models$models[[j]]
    cf <- coef(tmp)
    cf[] <- prm
    coef(tmp) <- cf

    ### check for fixed, offset, and shift_scale; go through predict if this is the case
    ### (slow but works)
    if (is.null(newdata)) {
        if (models$nn[j]) newdata <- tmp$data
    }

    if (models$cont[j]) {
        if (is.null(newdata)) {
            tr <- c(models$mm[[j]]$eY$Y %*% prm)
            trp <- c(models$mm[[j]]$eY$Yprime %*% prm)
            if (!models$normal[j])
                trd <- tmp$todistr$d(tr) * trp
        } else {
            tr <- predict(tmp, newdata = newdata, type = "trafo", ...)
            drv <- 1L
            names(drv) <- tmp$model$response
            trp <- predict(tmp, newdata = newdata, type = "trafo", 
                           deriv = drv, ...)
            if (!models$normal[j])
                trd <- predict(tmp, newdata = newdata, type = "density", ...)
        }
    } else {
        if (is.null(newdata)) {
            trl <- c(models$mm[[j]]$iY$Yleft %*% prm)
            trl[!is.finite(trl)] <- -Inf
            trr <- c(models$mm[[j]]$iY$Yright %*% prm)
            trr[!is.finite(trr)] <- Inf
        } else {
            stop("not yet implemented")
        }
    }

    if (what == "trafo") {
        stopifnot(models$cont[j])
        return(tr)
    }
    if (what == "dtrafo") {
        stopifnot(models$cont[j])
        return(tmp$todistr$d(tr))
    }
    if (what == "z") {
        stopifnot(models$cont[j])
        if (models$normal[j]) 
            return(tr)
        return(qnorm(tmp$todistr$p(tr, log = TRUE), log.p = TRUE))
    }
    if (what == "zleft") {
        stopifnot(!models$cont[j])
        if (models$normal[[j]])
            return(trl)
        return(qnorm(tmp$todistr$p(trl, log = TRUE), log.p = TRUE))
    }
    if (what == "dzleft") {
        stopifnot(!models$cont[j])
        if (models$normal[[j]])
            return(trl)
        return(tmp$todistr$d(trl))
    }
   if (what == "zright") {
        stopifnot(!models$cont[j])
        if (models$normal[[j]])
            return(trr)
        return(qnorm(tmp$todistr$p(trr, log = TRUE), log.p = TRUE))
    }
    if (what == "dzright") {
        stopifnot(!models$cont[j])
        if (models$normal[[j]])
            return(trr)
        return(tmp$todistr$d(trr))
    }
    if (what == "zprime") {
        stopifnot(models$cont[j])
        if (models$normal[[j]])
            return(trp)
        qn <- qnorm(tmp$todistr$p(tr, log = TRUE), log.p = TRUE)
        return(trd / dnorm(qn))
    }
    if (what == "estfun") {
        if (is.null(newdata))
            return(estfun(tmp))
        return(estfun(tmp, newdata = newdata))
    }
}

.MCw <- function(J, M, seed) {

    ### from stats:::simulate.lm
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
        runif(1)
    if (is.null(seed)) 
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }

    return(matrix(runif((J - 1) * M), ncol = M))
}

mmlt <- function(..., formula = ~ 1, data, conditional = FALSE, 
                 theta = NULL,
                 optim = mltoptim(auglag = list(maxtry = 5)), args = list(seed = 1, M = 1000), dofit = TRUE)
{
  
    call <- match.call()

    m <- .models(...)

    if (conditional && !all(m$normal))
        stop("Conditional models only available for marginal probit-type models")

    cJ <- sum(m$cont)
    dJ <- sum(!m$cont)
    J <- cJ + dJ
    Jp <- J * (J - 1) / 2
    llsc <- .ll(c(cJ, dJ), scale = !conditional, args)

    if (dJ && is.null(args$w))
        args$w <- .MCw(J = dJ, M = args$M, seed = args$seed)

    if (isTRUE(all.equal(formula, ~ 1))) {
        lX <- matrix(1)
        colnames(lX) <- "(Intercept)"
        bx <- NULL
    } else {
        bx <- formula
        if (inherits(formula, "formula")) {
            bx <- as.basis(formula, data)
        } 
        lX <- model.matrix(bx, data = data)
        if (conditional)
            warning("Conditional models with covariate-dependent correlations are order-dependent")
    }
    .Xparm <- function(parm) {
        parm <- parm[-(1:sum(m$nparm))]
        return(matrix(parm, nrow = ncol(lX)))
    }

    if (cJ) {
        mm <- lapply(1:cJ, function(j) .model_matrix(m, j = j))
        mmp <- lapply(1:cJ, function(j) .model_matrix(m, j = j, prime = TRUE))
    }
    if (dJ)
        dmm <- lapply(cJ + 1:dJ, function(j) .model_matrix(m, j = j))

    ### note: estfun() already has weights already multiplied to scores
    weights <- m$weights
    if (weights) {
        if (cJ) {
            mm <- lapply(mm, function(x) x * weights)
            mmp <- lapply(mmp, function(x) x * weights)
        }
        if (dJ) {
            dmm <- lapply(dmm, function(x) list(Yleft = x$Yleft * weights,
                                                Yright = x$Yright * weights))
        }
    }

    LAMBDA <- ltMatrices(matrix(0, nrow = Jp, ncol = nrow(lX)),
                         byrow = TRUE, diag = FALSE, names = names(m$models))

    ll <- function(parm) {

        # Lambda <- ltMatrices(t(lX %*% .Xparm(parm)), byrow = TRUE, diag = FALSE, 
        #                      names = names(m$models))
        # saves time in ltMatrices
        Lambda <- LAMBDA
        Lambda[] <- t(lX %*% .Xparm(parm))
        ret <- 0
        if (cJ) {
            z <- do.call("rbind", .mget(m, j = which(m$cont), parm = parm, what = "z"))
            zp <- do.call("rbind", .mget(m, j = which(m$cont), parm = parm, what = "zprime"))
            ret <- colSums(.log(zp))
            if (!dJ) return(ret + llsc$logLik(obs = z, Lambda = Lambda))
        }
        if (dJ) {
            lower <- do.call("rbind", .mget(m, j = which(!m$cont), parm = parm, what = "zleft"))
            upper <- do.call("rbind", .mget(m, j = which(!m$cont), parm = parm, what = "zright"))
            if (!cJ)
                return(llsc$logLik(lower = lower, upper = upper, Lambda = Lambda))
        }
        return(ret + llsc$logLik(obs = z, lower = lower, upper = upper, Lambda = Lambda))
    }

    sc <- function(parm) {

        # Lambda <- ltMatrices(t(lX %*% .Xparm(parm)), byrow = TRUE, diag = FALSE, 
        #                      names = names(m$models))
        # saves time in ltMatrices
        Lambda <- LAMBDA
        Lambda[] <- t(lX %*% .Xparm(parm))

        if (cJ) {
            z <- do.call("rbind", .mget(m, j = which(m$cont), parm = parm, what = "z"))
            if (!dJ)
                sc <- llsc$score(obs = z, Lambda = Lambda)
        }
        if (dJ) {
            lower <- do.call("rbind", .mget(m, j = which(!m$cont), parm = parm, what = "zleft"))
            upper <- do.call("rbind", .mget(m, j = which(!m$cont), parm = parm, what = "zright"))
            if (!cJ)
                sc <- llsc$score(lower = lower, upper = upper, Lambda = Lambda)
        }
        if (cJ && dJ)
            sc <- llsc$score(obs = z, lower = lower, upper = upper, Lambda = Lambda)

        Lmat <- Lower_tri(sc$Lambda)[rep(1:Jp, each = ncol(lX)), , drop = FALSE]
        if (identical(c(lX), 1)) {
            scL <- rowSums(Lmat)
        } else {
            scL <- rowSums(Lmat * t(lX[,rep(1:ncol(lX), Jp), drop = FALSE]))
        }
      
        scp <- vector(mode = "list", length = cJ + dJ)

        if (cJ) {
            if (all(m$normal)) {
                zp <- do.call("rbind", .mget(m, j = which(m$cont), parm = parm, what = "zprime"))
                scp[1:cJ] <- lapply(1:cJ, function(j) {
                    colSums(mm[[j]] * c(sc$obs[j,])) + colSums(mmp[[j]] / c(zp[j,]))
                })
            } else {
                dz <- do.call("rbind", .mget(m, j = which(m$cont), parm = parm, what = "dtrafo"))
                ef <- lapply(which(m$cont), function(j) .mget(m, j = j, parm = parm, what = "estfun"))
                scp[1:cJ] <- lapply(1:cJ, function(j) {
                    colSums(mm[[j]] * c(sc$obs[j,] + z[j,]) / c(dnorm(z[j,])) * c(dz[j,])) - colSums(ef[[j]])
                })
            }
        }

        if (dJ) {
            if (all(m$normal)) {
                scp[cJ + 1:dJ] <- lapply(1:dJ, function(j) {
                    colSums(dmm[[j]]$Yleft * c(sc$lower[j,])) +
                    colSums(dmm[[j]]$Yright * c(sc$upper[j,]))
                })
            } else {
                dzl <- do.call("rbind", .mget(m, j = which(!m$cont), parm = parm, what = "dzleft"))
                dzr <- do.call("rbind", .mget(m, j = which(!m$cont), parm = parm, what = "dzright"))
                scp[cJ + 1:dJ] <- lapply(1:dJ, function(j) {
                    dl <- c(dnorm(lower[j,]))
                    dl[!is.finite(lower[j,])] <- 1
                    dr <- c(dnorm(upper[j,]))
                    dr[!is.finite(upper[j,])] <- 1
                    colSums(dmm[[j]]$Yleft / dl * c(dzl[j,]) * c(sc$lower[j,])) +
                    colSums(dmm[[j]]$Yright / dr * c(dzr[j,]) * c(sc$upper[j,]))
                })
            }
        }
         
        ret <- c(do.call("c", scp), scL)
        return(ret)
    }

    if (is.null(theta)) {

        start <- do.call("c", lapply(m$models, function(mod) coef(mod)))
        if (weights) {
            cll <- function(cpar) -sum(weights * ll(c(start, cpar)))
        } else {
            cll <- function(cpar) -sum(ll(c(start, cpar)))
        }
        csc <- function(cpar) -sc(c(start, cpar))[-(1:length(start))]

        ### note: this is not optimal for conditional = TRUE
        for (i in 1:length(optim)) {
            op <- optim[[i]](theta = rep(0, Jp * ncol(lX)), f = cll, g = csc)
            if (op$convergence == 0) break()
        }
        # if (ret$convergence != 0)
        #     warning("Optimisation did not converge")

        start <- c(start, op$par)
    } else {
        ### use user-supplied starting values
        start <- theta
    }

    if (weights) {
        f <- function(par) -sum(weights * ll(par))
    } else {
        f <- function(par) -sum(ll(par))
    }
    g <- function(par) -sc(par)

    if (!dofit) return(list(ll = f, gr = g))

    ui <- m$ui
    ui <- cbind(ui, matrix(0, nrow = nrow(ui), ncol = Jp * ncol(lX)))
    ci <- m$ci
  
    for (i in 1:length(optim)) {
        ret <- optim[[i]](theta = start, f = f, g = g, ui = ui, ci = ci)
        if (ret$convergence == 0) break()
    }
    if (ret$convergence != 0)
        warning("Optimisation did not converge")

    pnm <- m$parm(ret$par)
    pnm <- sapply(1:J, function(j) paste(m$names[j], names(pnm[[j]]), sep = "."))
    tmp <- .Xparm(ret$par)
    rownames(tmp) <- colnames(lX)
    tmp <- unclass(ltMatrices(t(tmp), byrow = TRUE, diag = FALSE, names = m$names))
    tmp <- do.call("paste", expand.grid(rownames(tmp), colnames(tmp), sep = ".", 
                   stringsAsFactors = FALSE))
    names(ret$par) <- c(pnm, tmp)
  
    ret$ll <- f
    ret$sc <- g
    ret$args <- args
    ret$logLik <- -ret$value
    ret$models <- m
    ret$formula <- formula
    ret$bx <- bx
    ret$parm <- function(par) c(m$parm(par), list(.Xparm(par)))
    if (!missing(data))
        ret$data <- data
    ret$names <- m$names
    ret$call <- match.call()
    class(ret) <- c(ifelse(conditional, "cmmlt", "mmmlt"), "mmlt")
    ret
}


.coef.mmlt <- function(object, newdata, 
                      type = c("all", "Lambda", "Lambdainv", "Precision", "Sigma", "Corr", "Spearman"), 
                      ...)
{
  
    type <- match.arg(type)
    if (type == "all") return(object$par)

    if (type == "Spearman")
        return(6 * asin(coef(object, newdata = newdata, type = "Cor") / 2) / pi)
  
    prm <- object$parm(object$par)
    prm <- prm[[length(prm)]]

    if (missing(newdata) || is.null(object$bx)) {
        if (NROW(prm) > 1L && type != "Lambda")
            stop("newdata not specified")
        ret <- ltMatrices(t(prm), byrow = TRUE, diag = FALSE, names = object$names)
    } else {
        X <- model.matrix(object$bx, data = newdata)
        ret <- ltMatrices(t(X %*% prm), byrow = TRUE, diag = FALSE, names = object$names)
    }

    if (inherits(object, "mmmlt")) ret <- invcholD(ret)

    ret <- switch(type, "Lambda" = ret,
                        "Lambdainv" = solve(ret),
                        "Precision" = invchol2pre(ret),
                        "Sigma" = invchol2cov(ret),
                        "Corr" = invchol2cor(ret))
    return(ret)
}

coef.cmmlt <- function(object, newdata,
                       type = c("all", "conditional", "Lambda", "Lambdainv", 
                                "Precision", "Sigma", "Corr", "Spearman"), 
                       ...)
{

    type <- match.arg(type)
    if (type == "conditional") {
        prm <- object$parm(object$par)
        return(prm[-length(prm)])
    }
    return(.coef.mmlt(object = object, newdata = newdata, type = type, ...))
}

coef.mmmlt <- function(object, newdata,
                       type = c("all", "marginal", "Lambda", "Lambdainv", 
                                "Precision", "Sigma", "Corr", "Spearman"), 
                       ...)
{

    type <- match.arg(type)
    if (type == "marginal") {
        prm <- object$parm(object$par)
        return(prm[-length(prm)])
    }
    return(.coef.mmlt(object = object, newdata = newdata, type = type, ...))
}


vcov.mmlt <- function(object, ...) {
    step <- 0
    lam <- 1e-6
    H <- object$optim_hessian
    if (is.null(H)) {
        if (requireNamespace("numDeriv")) {
            H <- numDeriv::hessian(object$ll, object$par)
        } else {
            stop("Hessian not available")
        }
    }
    while((step <- step + 1) <= 3) {
          ret <- try(solve(H + (step - 1) * lam * diag(nrow(H))))
          if (!inherits(ret, "try-error")) break
    }
    if (inherits(ret, "try-error"))
        stop("Hessian is not invertible")
    if (step > 1)
        warning("Hessian is not invertible, an approximation is used")
    rownames(ret) <- colnames(ret) <- names(coef(object))
    ret
}

logLik.mmlt <- function (object, parm = coef(object), ...) 
{
    args <- list(...)
    if (length(args) > 0) 
        warning("Arguments ", names(args), " are ignored")
    ret <- -object$ll(parm)
    attr(ret, "df") <- length(object$par)
    class(ret) <- "logLik"
    ret
}

summary.mmlt <- function(object, ...) {
    ret <- list(call = object$call,
                #                tram = object$tram,
                test = cftest(object, parm = names(coef(object, with_baseline = FALSE))),
                ll = logLik(object))
    class(ret) <- "summary.mmlt"
    ret
}

print.summary.mmlt <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
    cat("\n", "Multivariate conditional transformation model", "\n")
    cat("\nCall:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    pq <- x$test$test
    mtests <- cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
    colnames(mtests) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    sig <- .Machine$double.eps
    printCoefmat(mtests, digits = digits, has.Pvalue = TRUE, 
                 P.values = TRUE, eps.Pvalue = sig)
    cat("\nLog-Likelihood:\n ", x$ll, " (df = ", attr(x$ll, "df"), 
          ")", sep = "")
    cat("\n\n")
    invisible(x)
}

print.mmlt <- function(x, ...) {
    cat("\n", "Multivariate conditional transformation model", "\n")
    cat("\nCall:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    print(coef(x))
    invisible(x)
}


predict.mmlt <- function (object, newdata, margins = 1:J, type = c("trafo", "distribution", 
    "density"), log = FALSE, ...) 
{
    type <- match.arg(type)
    J <- length(object$models$models)
    margins <- sort(margins)
    stopifnot(all(margins %in% 1:J))

    if (length(margins) == 1L) {
        ### ... may carry q = something
        tmp <- object$models$models[[margins]]
        cf <- coef(tmp)
        cf[] <- object$models$parm(coef(object))[[margins]]
        coef(tmp) <- cf
        ret <- predict(tmp, newdata = newdata, type = type, log = log, ...)
        return(ret)
    }

    ### don't feed ...
    z <- .mget(object$models, margins, parm = coef(object, type = "all"),
               newdata = newdata, what = "z")
    z <- do.call("rbind", z)

    if (type == "trafo") {
        stopifnot(!log)
        L <- coef(object, newdata = newdata, type = "Lambda")
        if (length(margins) != J) 
            L <- marg_mvnorm(invchol = L, which = margins)$invchol
        return(Mult(L, z))
    }
    if (type == "distribution") {
        upper <- z
        lower <- matrix(-Inf, ncol = ncol(z), nrow = nrow(z))
        Linv <- coef(object, newdata = newdata, type = "Lambdainv")
        if (length(margins) != J) 
            Linv <- marg_mvnorm(chol = Linv, which = margins)$chol
        a <- object$args
        a$lower <- lower
        a$upper <- upper
        a$logLik <- FALSE
        a$chol <- Linv
        ret <- do.call("lpmvnorm", a)
        if (log) return(ret)
        return(exp(ret))
    }
    stopifnot(all(object$models$cont))
    zprime <- .mget(object$models, margins, parm = coef(object, type = "all"),
                    newdata = newdata, what = "zprime")
    if (length(margins) > 1L) {
        zprime <- do.call("rbind", zprime)
    } else {
        zprime <- matrix(zprime, nrow = 1)
    }
    L <- coef(object, newdata = newdata, type = "Lambda")
        if (length(margins) != J) 
            L <- marg_mvnorm(invchol = L, which = margins)$invchol
    ret <- dmvnorm(t(z), invchol = L, log = TRUE)
    ret <- ret + colSums(.log(zprime))
    if (log) return(ret)
    return(exp(ret))
}

simulate.mmlt <- function(object, nsim = 1L, seed = NULL, newdata, K = 50, ...) {

    ### from stats:::simulate.lm
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
        runif(1)
    if (is.null(seed)) 
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }

    if (!is.data.frame(newdata))
        stop("not yet implemented")

    args <- list(...)
    if (length(args) > 0L)
        stop("argument(s)", paste(names(args), collapse = ", "), "ignored")

    if (nsim > 1L) 
        return(replicate(nsim, simulate(object, newdata = newdata, K = K, ...), 
                         simplify = FALSE))

    J <- length(object$models$models)
    L <- coef(object, newdata = newdata, type = "Lambda")
    N <- nrow(newdata)

    Z <- matrix(rnorm(J * N), ncol = N)
    Ztilde <- solve(L, Z)

    ret <- matrix(0.0, nrow = N, ncol = J)

    if (inherits(object, "cmmlt")) {
        for (j in 1:J) {
            q <- mkgrid(object$models$models[[j]], n = K)[[1L]]
            pr <- predict(object$models$models[[j]], newdata = newdata, type = "trafo", q = q)
            if (!is.matrix(pr)) pr <- matrix(pr, nrow = length(pr), ncol = NROW(newdata))
            ret[,j] <- as.double(mlt:::.invf(object$models$models[[j]], f = t(pr), 
                                             q = q, z = t(Ztilde[j,,drop = FALSE])))
        }
    } else {
        Ztilde <- pnorm(Ztilde, log.p = TRUE)
        for (j in 1:J) {
            q <- mkgrid(object$models$models[[j]], n = K)[[1L]]
            pr <- predict(object$models$models[[j]], newdata = newdata, type = "logdistribution", q = q)
            if (!is.matrix(pr)) pr <- matrix(pr, nrow = length(pr), ncol = NROW(newdata))
            ret[,j] <- as.double(mlt:::.invf(object$models$models[[j]], f = t(pr), 
                                             q = q, z = t(Ztilde[j,,drop = FALSE])))
        }
    }
    colnames(ret) <- variable.names(object, response_only = TRUE)
    return(ret)
}

variable.names.mmlt <- function(object, response_only = FALSE, ...) {

    if (response_only)
        return(sapply(object$models$models, function(x) variable.names(x)[1L]))
    vn <- unique(c(sapply(object$models$models, function(x) variable.names(x)), 
                 all.vars(object$formula)))
    return(vn)
}
    
confregion <- function(object, level = .95, ...)
    UseMethod("confregion")

confregion.mmlt <- function(object, level = .95, newdata, K = 250, ...) {

    if (!missing(newdata)) stopifnot(nrow(newdata) == 1)

    if (object$conditional) {
        Linv <- coef(object, newdata = newdata, type = "Lambdainv")
        Linv <- as.array(Linv)[,,1]
    } else {
        CR <- as.array(coef(object, newdata = newdata, type = "Corr"))[,,1]
        Linv <- t(chol(CR))
    }
    J <- nrow(Linv)

    q <- qchisq(level, df = J)

    if (J == 2) {
        angle <- seq(0, 2 * pi, length = K)
        x <- cbind(cos(angle), sin(angle))
    } else {
        x <- matrix(rnorm(K * J), nrow = K, ncol = J)
        x <- x / sqrt(rowSums(x^2))
    }
    x <- sqrt(q) * x
    a <- x %*% t(Linv)

    nd <- if (missing(newdata)) data.frame(1) else newdata

    ret <- lapply(1:J, function(j) {
        prb <- object$models$models[[j]]$todistr$p(a[,j])
        predict(object$models$models[[j]], newdata = nd, type = "quantile", prob = prb)
    })
    
    ret <- do.call("cbind", ret)
    return(ret)
}

HDR <- function(object, level = .95, ...)
    UseMethod("HDR")

HDR.mmlt <- function(object, level = .95, newdata, nsim = 1000L, K = 25, ...) {

    if (!missing(newdata)) {
        stopifnot(nrow(newdata) == 1)
    } else {
        newdata <- data.frame(1)
    }

    ### https://doi.org/10.2307/2684423 Section 3.2
    y <- simulate(object, newdata = newdata[rep(1, nsim),,drop = FALSE])
    y <- cbind(y, newdata)
    d <- predict(object, newdata = y, type = "density")

    ret <- do.call("expand.grid", lapply(object$models$models, function(x) mkgrid(x, n = K)[[1L]]))
    colnames(ret) <- variable.names(object, response_only = TRUE)
    ret <- cbind(ret, newdata)
    ret$density <- predict(object, newdata = ret, type = "density")
    attr(ret, "cuts") <- quantile(d, prob = 1 - level)
    ret
}

mkgrid.mmlt <- function(object, ...) {

    lx <- mkgrid(as.basis(object$formula, data = object$data), ...)
    grd <- do.call("c", lapply(object$models$models, mkgrid, ...))
    grd <- c(grd, lx)
    do.call("expand.grid", grd[unique(names(grd))])
}
