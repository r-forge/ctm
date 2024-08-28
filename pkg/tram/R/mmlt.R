
### catch constraint violations here
.log <- function(x) {
  ret <- log(pmax(.Machine$double.eps, x))
  dim(ret) <- dim(x)
  ret
}

.rbind <- function(x) {
    if (!is.list(x)) return(matrix(x, nrow = 1))
    return(do.call("rbind", x))
}

.ll <- function(dim, standardize = TRUE, args = list()) {

    if (length(dim) == 1L)
        dim <- c(dim, 0L)

    cJ <- dim[1L]
    dJ <- dim[2L]

    if (!dJ) {

        cll <- function(obs, Lambda) {

            if (dim(Lambda)[2L] > 1)
                stopifnot(!attr(Lambda, "diag"))

            return(logLik(mvnorm(invchol = Lambda), obs = obs, 
                          standardize = standardize, logLik = FALSE))
        }

        csc <- function(obs, Lambda) {

            if (dim(Lambda)[2L] > 1)
                stopifnot(!attr(Lambda, "diag"))

            ret <- lLgrad(mvnorm(invchol = Lambda), obs = obs, standardize = standardize)
            return(list(Lambda = ret$scale, obs = ret$obs))
       }

       return(list(logLik = cll, score = csc))
    }

    ll <- function(obs = NULL, lower, upper, Lambda) {

        if (dim(Lambda)[2L] > 1)
            stopifnot(!attr(Lambda, "diag"))

        a <- args
        a$object <- mvnorm(invchol = Lambda)
        a$obs <- obs
        a$lower <- lower
        a$upper <- upper
        a$standardize <- standardize
        a$logLik <- FALSE
        return(do.call("logLik", a))
    }

    sc <- function(obs = NULL, lower, upper, Lambda) {

        a <- args
        a$object <- mvnorm(invchol = Lambda)
        a$obs <- obs
        a$lower <- lower
        a$upper <- upper
        a$standardize <- standardize
        ret <- do.call("lLgrad", a)
        ret <- list(Lambda = ret$scale,
                    obs = ret$obs,
                    mean = ret$mean, 
                    lower = ret$lower, 
                    upper = ret$upper)
        return(ret)
    }

    return(list(logLik = ll, score = sc))
}

.models <- function(...) {

    m <- lapply(list(...), function(x) as.mlt(x))
    # nm <- abbreviate(sapply(m, function(x) x$model$response), 4)
    nm <- sapply(m, function(x) x$model$response)
    J <- length(m)
    Jp <- J * (J - 1) / 2
    normal <- sapply(m, function(x) x$todistr$name == "normal")
  
    w <- lapply(m, weights)
    out <- lapply(w, function(x) stopifnot(isTRUE(all.equal(x, w[[1]]))))
    w <- w[[1L]]
  
    ### determine if response is numeric and was measured exactly
    ### (censoring and discreteness are treated the same here)
    mm <- lapply(m, function(mod) {
      eY <- get("eY", environment(mod$parm))
      iY <- get("iY", environment(mod$parm))
      list(eY = eY, iY = iY)
    })

    cmod <- sapply(mm, function(x) !is.null(x$eY))  
    dmod <- sapply(mm, function(x) !is.null(x$iY))  
    stopifnot(all(xor(cmod, dmod)))
    ### <FIXME> is this the only place (except in ldpmvnorm) 
    ### where this is checked and assumed?
    ### </FIXME>
    ### continuous models first
#    stopifnot(all(diff(cmod) <= 0))
#    stopifnot(all(diff(dmod) >= 0))

    ### determine if response is conceptually numeric
    cresp <- sapply(m, function(x) 
        inherits(attr(x$model$bases$response, "variables"), 
                 "continuous_var"))

    nobs <- unique(sapply(m, nobs))
    stopifnot(length(nobs) == 1L)
    nobs <- nobs[[1L]]

    P <- sapply(m, function(x) length(coef(x, fixed = TRUE)))
    fpar <- factor(rep(1:J, P))

    ### par always includes marginally fixed parameters
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

    fixed <- vector(mode = "list", length = J)
    for (j in 1:J) {
        if (!is.null(m[[j]]$fixed)) {
            fj <- m[[j]]$fixed
            names(fj) <- paste(nm[j], names(fj), sep = ".")
            fixed[[j]] <- fj
        }
    }

    return(list(models = m, mf = mf, cont = cmod, cresp = cresp, 
                normal = normal, nobs = nobs, weights = w, 
                nparm = P, parm = parm, ui = ui, ci = ci, mm = mm, 
                names = nm, fixed = fixed))
}

.mget <- function(models, j = 1, parm, newdata = NULL, weights = NULL,
                  what = c("trafo", "dtrafo", "z", "zleft", 
                           "dzleft", "zright", "dzright", "zprime", 
                           "trafoprime", "estfun", "scale"), ...) {

    what <- match.arg(what)

    if (length(j) > 1) {
        ret <- lapply(j, .mget, models = models, parm = parm, 
                      newdata = newdata, weights = weights, 
                      what = what, ...)
        return(ret)
    }

    if (what == "scale")
        return(models$models[[j]]$parsc)

    prm <- models$parm(parm)[[j]]
    ### remove marginally fix parameters
    if (!is.null(models$fixed[[j]]))
        prm <- prm[!names(prm) %in% names(models$fixed[[j]])]
    tmp <- as.mlt(models$models[[j]])
    if (!is.null(newdata)) {
        tmp <- mlt(tmp$model, data = newdata, # weights = weights,
                                              # offset = tmp$offset, 
                   fixed = tmp$fixed, theta = prm,
                   scale = tmp$scale, dofit = FALSE)
    }


    if (what == "trafoprime") {
        ret <- models$models[[j]]$trafoprime(prm, weights)
        if (models$cont[j])
            return(ret[c("exY", "exYprime")])
        ret$iYleft[!is.finite(ret$iYleft[,1])] <- 0
        ret$iYright[!is.finite(ret$iYright[,1])] <- 0
        return(ret[c("iYleft", "iYright")])
    }

    ret <- tmp$trafo(prm, weights)
    ### extract both exact and interval (former might be needed for
    ### predictions)
    tr <- ret$trex
    trp <- ret$trexprime
    if (!models$normal[j])
        trd <- tmp$todistr$d(tr) * trp
    trl <- ret$trleft
    trr <- ret$trright

    if (what == "trafo") {
        return(tr)
    }
    if (what == "dtrafo") {
        return(tmp$todistr$d(tr))
    }
    if (what == "z") {
        if (models$normal[j]) 
            return(tr)
        return(qnorm(tmp$todistr$p(tr, log = TRUE), log.p = TRUE))
    }
    if (what == "zleft") {
        if (models$normal[[j]])
            return(trl)
        return(qnorm(tmp$todistr$p(trl, log = TRUE), log.p = TRUE))
    }
    if (what == "dzleft") {
        if (models$normal[[j]])
            return(rep(1, length(trl)))
        qn <- qnorm(tmp$todistr$p(trl, log = TRUE), log.p = TRUE)
        dn <- dnorm(qn)
        dn[!is.finite(dn)] <- 1
        return(tmp$todistr$d(trl) / dn)
    }
   if (what == "zright") {
        if (models$normal[[j]])
            return(trr)
        return(qnorm(tmp$todistr$p(trr, log = TRUE), log.p = TRUE))
    }
    if (what == "dzright") {
        if (models$normal[[j]])
            return(rep(1, length(trr)))
        qn <- qnorm(tmp$todistr$p(trr, log = TRUE), log.p = TRUE)
        dn <- dnorm(qn)
        dn[!is.finite(dn)] <- 1
        return(tmp$todistr$d(trr) / dn)
    }
    if (what == "zprime") {
        if (models$normal[[j]])
            return(trp)
        qn <- qnorm(tmp$todistr$p(tr, log = TRUE), log.p = TRUE)
        return(trd / dnorm(qn))
    }
    if (what == "estfun") {
        return(estfun(tmp, parm = prm))
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

.start <- function(m, xnames, fixed = NULL) {

    J <- length(m$models)
    Jp <- J * (J - 1) / 2
    Jnames <- m$names
    ### fixed = FALSE?
    margin_par <- do.call("c", lapply(m$models, 
                                      function(mod) coef(as.mlt(mod))))
    names(margin_par) <- paste(rep(Jnames, time = m$nparm), 
                               names(margin_par), sep = ".")

    rn <- rownames(unclass(ltMatrices(1:Jp, names = Jnames, byrow = TRUE)))
    lnames <- paste(rep(rn, each = length(xnames)),
                    rep(xnames, length(rn)), sep = ".")
    lambda_par <- rep(0, length(lnames))
    names(lambda_par) <- lnames

    start <- c(margin_par, lambda_par)

    if (!is.null(fixed))
        stopifnot(all(fixed %in% names(start)))

    return(start)
}

mmltoptim <- function(auglag = list(maxtry = 5), ...)
    mltoptim(auglag = auglag, ...)

mmlt <- function(..., formula = ~ 1, data, conditional = FALSE, 
                 theta = NULL, fixed = NULL, scale = FALSE,
                 optim = mmltoptim(), args = list(seed = 1, M = 1000), 
                 dofit = TRUE, domargins = TRUE, sequentialfit = FALSE)
{
  
    call <- match.call()

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
            warning("Conditional models with covariate-dependent",
                    "correlations are order-dependent")
    }

    if (sequentialfit) {
        m <- list(...)
        if (!dofit)
            stop("Cannot perform sequential fit")
        if (!is.null(theta))
            stop("Cannot perform sequential fit with starting values")
        if (!is.null(fixed))
            stop("Cannot perform sequential fit with fixed values")
        mj <- as.mlt(m[[1]])
        for (j in 2:length(m)) {
            mc <- m[1:j]
            mc$formula <- formula
            mc$data <- data
            mc$conditional <- conditional
            mc$scale <- scale
            mc$optim <- optim
            mc$args <- args
            if (!is.null(mc$args$w)) 
                mc$args$w <- mc$args$w[1:(j - 1L),,drop = FALSE]
            mc$domargins <- domargins
            cf <- .start(do.call(".models", m[1:j]), colnames(lX))
            cfj <- coef(mj, fixed = TRUE)
            if (j == 2)
                names(cfj) <- names(cf)[1:length(cfj)]
            ### fix coefs and estimate jth row of
            ### lambda and (if domargins) jth marginal parameters only
            mc$fixed <- cfj
            mj <- do.call("mmlt", mc)
        }
        return(mj)
    }

    m <- .models(...)

    if (conditional && !all(m$normal))
        stop("Conditional models only available", 
             "for marginal probit-type models.")

    if (conditional && !domargins)
        stop("Conditional models must fit marginal and joint parameters.")

    ### compute starting values for lambda
    if (is.null(theta) && dofit && domargins) {
        cl <- match.call()
        cl$conditional <- FALSE
        cl$domargins <- FALSE
        sm <- eval(cl, parent.frame())
        if (!is.null(sm)) {
            theta <- coef(sm, type = "all", fixed = TRUE)
            if (conditional) {
                ### theta are conditional parameters, scale with sigma
                class(sm)[1] <- "cmmlt" ### do NOT standardize Lambda
                d <- rowMeans(mvtnorm::diagonals(coef(sm, newdata = data, 
                                                 type = "Sigma")))
                theta[1:sum(m$nparm)] <- 
                    theta[1:sum(m$nparm)] * rep(sqrt(d), times = m$nparm)
            }
        } else {
            theta <- do.call("c", lapply(m$models, function(mod) coef(mod)))
        }
    } 

    cJ <- sum(m$cont)
    dJ <- sum(!m$cont)
    J <- cJ + dJ
    Jp <- J * (J - 1) / 2
    llsc <- .ll(c(cJ, dJ), standardize = !conditional, args)

    if (dJ > 1L) {
        if (is.null(args$w)) {
            args$w <- .MCw(J = dJ, M = args$M, seed = args$seed)
        } else {
            if (!is.matrix(args$w)) args$w <- matrix(args$w, nrow = 1)
            if (nrow(args$w) < dJ - 1) stop("incorrect dimension of w")
            if (nrow(args$w) > dJ - 1)
                ### make sure only dJ - 1 columns are present
                args$w <- args$w[-(dJ:nrow(args$w)),,drop = FALSE]
        }
    }

    .Xparm <- function(parm) {
        parm <- parm[-(1:sum(m$nparm))]
        return(matrix(parm, nrow = ncol(lX)))
    }

    start <- .start(m, colnames(lX), names(fixed))
    ### marginal fixed parameters
    if (!is.null(m$fixed)) {
        mfixed <- do.call("c", m$fixed)
        fixed <- c(fixed, mfixed)
        start <- .start(m, colnames(lX), names(fixed))
    }
    
    parnames <- eparnames <- names(start)
    lparnames <- names(start)[-(1:sum(m$nparm))]
    if (!is.null(fixed)) eparnames <- eparnames[!eparnames %in% names(fixed)]
    if (!is.null(fixed)) lparnames <- lparnames[!lparnames %in% names(fixed)]

    if (!is.null(theta)) {
        if (!is.null(fixed)) theta <- theta[!names(theta) %in% names(fixed)]
        stopifnot(length(theta) == length(eparnames))
        names(theta) <- eparnames
    }

    ### note: estfun() already has weights already multiplied to scores
    weights <- m$weights

    LAMBDA <- ltMatrices(matrix(0, nrow = Jp, ncol = nrow(lX)),
                         byrow = TRUE, diag = FALSE, names = m$names) #names(m$models))

    ll <- function(parm, newdata = NULL) {

        if (!is.null(newdata) && !isTRUE(all.equal(formula, ~ 1))) 
            lX <- model.matrix(bx, data = newdata)

        # Lambda <- ltMatrices(t(lX %*% .Xparm(parm)), byrow = TRUE, 
        #                      diag = FALSE, names = m$names) ##names(m$models))
        # saves time in ltMatrices
        Lambda <- LAMBDA
        Lambda[] <- t(lX %*% .Xparm(parm))
        ret <- 0
        if (cJ) {
            z <- .rbind(.mget(m, j = which(m$cont), parm = parm, what = "z", 
                              newdata = newdata, weights = weights))
            rownames(z) <- m$names[which(m$cont)]
            zp <- .rbind(.mget(m, j = which(m$cont), parm = parm, 
                               what = "zprime", newdata = newdata, 
                               weights = weights))
            ret <- colSums(.log(zp))
            if (!dJ) return(ret + llsc$logLik(obs = z, Lambda = Lambda))
        }
        if (dJ) {
            lower <- .rbind(.mget(m, j = which(!m$cont), parm = parm, 
                                  what = "zleft", newdata = newdata, 
                                  weights = weights))
            upper <- .rbind(.mget(m, j = which(!m$cont), parm = parm, 
                                  what = "zright", newdata = newdata, 
                                  weights = weights))
            rownames(lower) <- rownames(upper) <- m$names[which(!m$cont)]
            if (!cJ)
                return(llsc$logLik(lower = lower, upper = upper, 
                                   Lambda = Lambda))
        }
        return(ret + llsc$logLik(obs = z, lower = lower, upper = upper, 
                           Lambda = Lambda))
    }

    sc <- function(parm, newdata = NULL, scores = FALSE) {

        if (scores) {
            RS <- CS <- function(x) x
        } else {
            RS <- function(x) rowSums(x, na.rm = TRUE)
            CS <- function(x) colSums(x, na.rm = TRUE)
        }

        if (!is.null(newdata) && !isTRUE(all.equal(formula, ~ 1))) 
            lX <- model.matrix(bx, data = newdata)

        # Lambda <- ltMatrices(t(lX %*% .Xparm(parm)), byrow = TRUE, 
        #                      diag = FALSE, names = m$names) # names(m$models))
        # saves time in ltMatrices
        Lambda <- LAMBDA
        Lambda[] <- t(lX %*% .Xparm(parm))

        if (cJ) {
            z <- .rbind(.mget(m, j = which(m$cont), parm = parm, what = "z", 
                              newdata = newdata, weights = weights))
            rownames(z) <- m$names[which(m$cont)]
            if (!dJ)
                sc <- llsc$score(obs = z, Lambda = Lambda)
        }
        if (dJ) {
            lower <- .rbind(.mget(m, j = which(!m$cont), parm = parm, 
                                  what = "zleft", newdata = newdata, 
                                  weights = weights))
            upper <- .rbind(.mget(m, j = which(!m$cont), parm = parm, 
                                  what = "zright", newdata = newdata, 
                                  weights = weights))
            rownames(lower) <- rownames(upper) <- m$names[which(!m$cont)]
            if (!cJ)
                sc <- llsc$score(lower = lower, upper = upper, 
                                 Lambda = Lambda)
        }
        if (cJ && dJ)
            sc <- llsc$score(obs = z, lower = lower, upper = upper, 
                             Lambda = Lambda)

        ### <FIXME> explain subset </FIXME>
        Lmat <- Lower_tri(sc$Lambda)[rep(1:Jp, each = ncol(lX)), , drop = FALSE]
        if (identical(c(lX), 1)) {
            scL <- RS(Lmat) ### NaN might appear in scores
        } else {
            scL <- RS(Lmat * t(lX[,rep(1:ncol(lX), Jp), drop = FALSE]))
        }
      
        scp <- vector(mode = "list", length = cJ + dJ)

        if (cJ) {
            mm <- lapply(which(m$cont), 
                function(j) .mget(m, j = j, parm = parm, what = "trafoprime", 
                                  newdata = newdata, weights = weights))
            if (all(m$normal)) {
                zp <- .rbind(.mget(m, j = which(m$cont), parm = parm, 
                                   what = "zprime", newdata = newdata, 
                                   weights = weights))
                scp[1:cJ] <- lapply(1:cJ, function(j) {
                    CS(mm[[j]]$exY * c(sc$obs[j,])) + 
                        CS(mm[[j]]$exYprime / c(zp[j,]))
                })
            } else {
                dz <- .rbind(.mget(m, j = which(m$cont), parm = parm, 
                                   what = "dtrafo", newdata = newdata, 
                                   weights = weights))
                ef <- lapply(which(m$cont), 
                             function(j) 
                                 .mget(m, j = j, parm = parm, what = "estfun", 
                                       newdata = newdata, weights = weights))
                scp[1:cJ] <- lapply(1:cJ, function(j) {
                    CS(mm[[j]]$exY * c(sc$obs[j,] + z[j,]) / 
                        c(dnorm(z[j,])) * c(dz[j,])) - CS(ef[[j]])
                })
            }
        }

        if (dJ) {
            mm <- lapply(which(!m$cont), 
                function(j) .mget(m, j = j, parm = parm, what = "trafoprime", 
                                  newdata = newdata, weights = weights))
            if (all(m$normal)) {
                scp[cJ + 1:dJ] <- lapply(1:dJ, function(j) {
                    CS(mm[[j]]$iYleft * c(sc$lower[j,])) +
                    CS(mm[[j]]$iYright * c(sc$upper[j,]))
                })
            } else {
                dzl <- .rbind(.mget(m, j = which(!m$cont), parm = parm, 
                                    what = "dzleft", newdata = newdata, 
                                    weights = weights))
                dzl[!is.finite(dzl)] <- 0
                dzr <- .rbind(.mget(m, j = which(!m$cont), parm = parm, 
                                    what = "dzright", newdata = newdata, 
                                    weights = weights))
                dzr[!is.finite(dzr)] <- 0
                scp[cJ + 1:dJ] <- lapply(1:dJ, function(j) {
                    return(CS(mm[[j]]$iYleft * c(dzl[j,]) * c(sc$lower[j,])) +
                           CS(mm[[j]]$iYright * c(dzr[j,]) * c(sc$upper[j,])))
                })
            }
        }
        
        if (!scores) {
            ret <- c(do.call("c", scp), c(scL))
            names(ret) <- parnames
            return(ret)
        }
        ret <- cbind(do.call("cbind", scp), t(scL))
        colnames(ret) <- parnames
        return(ret)
    }

    ### scale
    if (scale) {
        scl <- rep(apply(abs(lX), 2, max, na.rm = TRUE), times = Jp)
        lt1 <- scl < 1.1
        gt1 <- scl >= 1.1
        scl[gt1] <- 1 / scl[gt1]
        scl[lt1] <- 1
        scl <- c(do.call("c", .mget(m, j = 1:J, parm = NULL, what = "scale")), 
                              scl)
        names(scl) <- parnames
    } else {
        scl <- numeric(length(parnames))
        scl[] <- 1
        names(scl) <- parnames
    }

    f <- function(par, scl, ...) {
        if (!is.null(fixed)) {
            p <- par
            names(p) <- eparnames
            p <- c(p, fixed)
            par <- p[parnames]
        }
        return(-sum(weights * ll(par * scl, ...)))
    }

    ### note: x * weights was already computed
    g <- function(par, scl, ...) {
        if (!is.null(fixed)) {
            p <- par
            names(p) <- eparnames
            p <- c(p, fixed)
            par <- p[parnames]
        }
        ret <- -sc(par * scl, ...) * scl
        if (is.null(fixed)) return(ret)
        if (is.matrix(ret))
            return(ret[, !parnames %in% names(fixed)])
        return(ret[!parnames %in% names(fixed)])
    }

    if (!domargins) {

        stopifnot(!conditional)

        start <- start[1:sum(m$nparm)]
        start <- start[names(start) %in% eparnames]

        cll <- function(cpar) f(c(start / scl[names(start)], cpar), scl = scl)
        csc <- function(cpar) {
            ret <- g(c(start / scl[names(start)], cpar), scl = scl)
            return(ret[names(lambdastart)])
        }

        if (is.null(theta)) {
            lambdastart <- rep(0, length(lparnames))
            names(lambdastart) <- lparnames
        } else {
            lambdastart <- theta[lparnames]
        }

        if (length(lambdastart)) {
            for (i in 1:length(optim)) {
                op <- optim[[i]](theta = lambdastart, f = cll, g = csc)
                if (op$convergence == 0) break()
            }
            names(op$par) <- names(lambdastart)
            ret <- c(start, op$par * scl[names(lambdastart)])
            ### note: We throw away optim_hessian. vcov.mmlt
            ### uses numDeriv::hessian INCLUDING marginal parameters
            ### (which is probably the right thing to do)
            ret <- list(par = ret, value = -op$value)
        } else {
            ### no parameters to optimise over
            return(NULL)
        }
    } else {

        ui <- m$ui
        ui <- cbind(ui, matrix(0, nrow = nrow(ui), 
                               ncol = length(parnames) - ncol(ui)))
        ci <- m$ci
        if (!is.null(fixed)) {
            d <- ui[, parnames %in% names(fixed), drop = FALSE] %*% fixed
            ui <- ui[, !parnames %in% names(fixed), drop = FALSE]
            ci <- m$ci - d
        } 

        if (is.null(theta) && !dofit) 
            return(list(ll = function(...) f(..., scl = 1), 
                        score = function(...) g(..., scl = 1), 
                        ui = ui, ci = ci))

        start <- theta / scl[eparnames]
        ui <- t(t(ui) * scl[eparnames])
  
        if (dofit) {
            for (i in 1:length(optim)) {
                ret <- optim[[i]](theta = start, 
                                  f = function(par) f(par, scl = scl), 
                                  g = function(par) g(par, scl = scl), 
                                  ui = ui, ci = ci)
                if (ret$convergence == 0) break()
            }
            if (ret$convergence != 0)
                warning("Optimisation did not converge")
        } else {
            ret <- list(par = start, value = f(theta, scl = 1), 
                        convergence = NA, optim_hessian = NA)
        }
        names(ret$par) <- eparnames
        ret$par[eparnames] <- ret$par[eparnames] * scl[eparnames]
    }
  
    ret$ll <- function(...) f(..., scl = 1)
    ret$score <- function(...) g(..., scl = 1)
    ret$args <- args
    ret$logLik <- -ret$value
    ret$models <- m
    ret$formula <- formula
    ret$bx <- bx
    ret$parm <- function(par, flat = FALSE) {
        if (!is.null(fixed)) 
            par <- c(par, fixed)[parnames]
        if (flat) return(par)
        return(c(m$parm(par), list(.Xparm(par))))
    }
    if (!missing(data))
        ret$data <- data
    ret$names <- m$names
    ret$call <- match.call()
    class(ret) <- c(ifelse(conditional, "cmmlt", "mmmlt"), "mmlt")
    ret$mmlt <- "Multivariate Conditional Transformation Model"
    ret
}


.coef.mmlt <- function(object, newdata,
                       type = c("all", "Lambdapar", "Lambda", "Lambdainv", 
                                "Precision", "PartialCorr", "Sigma", "Corr", 
                                "Spearman", "Kendall"), 
                       fixed = FALSE, ...)
{
  
    type <- match.arg(type)
    if (type == "all") {
        if (!fixed) return(object$par)
        return(object$parm(object$par, flat = TRUE))
    }

    if (type == "Spearman")
        return(6 * asin(coef(object, newdata = newdata, type = "Cor") / 2) / pi)
  
    if (type == "Kendall")
        return(2 * asin(coef(object, newdata = newdata, type = "Cor")) / pi)

    prm <- object$parm(object$par)
    prm <- prm[[length(prm)]]

    if (missing(newdata) || is.null(object$bx)) {
        if (NROW(prm) > 1L && type != "Lambda")
            stop("newdata not specified")
        ret <- ltMatrices(t(prm), byrow = TRUE, diag = FALSE, 
                          names = object$names)
    } else {
        X <- model.matrix(object$bx, data = newdata)
        ret <- ltMatrices(t(X %*% prm), byrow = TRUE, diag = FALSE, 
                          names = object$names)
    }

    if (inherits(object, "mmmlt")) {
        ret0 <- ret
        ret <- mvtnorm::standardize(invchol = ret)
    }

    ret <- switch(type, "Lambdapar" = ret0, 
                        "Lambda" = ret,
                        "Lambdainv" = solve(ret),
                        "Precision" = invchol2pre(ret),
                        "PartialCorr" = invchol2pc(ret),
                        "Sigma" = invchol2cov(ret),
                        "Corr" = invchol2cor(ret))
    return(ret)
}

coef.cmmlt <- function(object, newdata,
                       type = c("all", "conditional", "Lambdapar", "Lambda", 
                                "Lambdainv", "Precision", "PartialCorr", 
                                "Sigma", "Corr", "Spearman", "Kendall"), 
                       fixed = FALSE,
                       ...)
{

    type <- match.arg(type)
    if (type == "conditional") {
        prm <- object$parm(object$par)
        return(prm[-length(prm)])
    }
    return(.coef.mmlt(object = object, newdata = newdata, type = type, 
                      fixed = fixed, ...))
}

coef.mmmlt <- function(object, newdata,
                       type = c("all", "marginal", "Lambdapar", "Lambda", 
                                "Lambdainv", "Precision", "PartialCorr", 
                                "Sigma", "Corr", "Spearman", "Kendall"), 
                       fixed = FALSE,
                       ...)
{

    type <- match.arg(type)
    if (type == "marginal") {
        prm <- object$parm(object$par)
        return(prm[-length(prm)])
    }
    return(.coef.mmlt(object = object, newdata = newdata, type = type, 
                      fixed = fixed, ...))
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

    ### <NOTE> add an option to compute vcov for selected 
    ### parameters (eg marginal effects) only and use Schur
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

logLik.mmlt <- function (object, parm = coef(object), newdata = NULL, ...) 
{
    args <- list(...)
    if (length(args) > 0) 
        warning("Arguments ", names(args), " are ignored")
    ret <- -object$ll(parm, newdata = newdata)
    attr(ret, "df") <- length(object$par)
    class(ret) <- "logLik"
    ret
}

estfun.mmlt <- function(x, parm = coef(x, type = "all"), 
                        newdata = NULL, ...) {
    args <- list(...)
    if (length(args) > 0)
        warning("Arguments ", names(args), " are ignored")
    return(x$score(parm, newdata = newdata, scores = TRUE))
}

summary.mmlt <- function(object, ...) {
    ret <- list(call = object$call,
                #                tram = object$tram,
                test = cftest(object, 
                              parm = names(coef(object, with_baseline = FALSE))),
                ll = logLik(object))
    class(ret) <- "summary.mmlt"
    ret
}

print.summary.mmlt <- function(x, digits = max(3L, getOption("digits") - 3L), 
                               ...) {
    cat("\n", x$mmlt, "\n")
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


predict.mmlt <- function (object, newdata, margins = 1:J, 
    type = c("trafo", "distribution", "survivor", "density", "hazard"), 
    log = FALSE, args = object$args, ...) 
{
    J <- length(object$models$models)
    margins <- sort(margins)
    stopifnot(all(margins %in% 1:J))

    if (length(margins) == 1L) {
        ### ... may carry q = something
        tmp <- object$models$models[[margins]]
        cf <- coef(tmp, fixed = TRUE)
        ncf <- names(cf)
        names(cf) <- paste(variable.names(tmp)[1L], names(cf), sep = ".")
        cfm <- object$models$parm(coef(object, fixed = TRUE))[[margins]]
        cf[names(cfm)] <- cfm
        names(cf) <- ncf
        coef(tmp) <- cf
        ### marginal models
        if (!inherits(object, "cmmlt")) {
            ret <- predict(tmp, newdata = newdata, type = type, log = log, ...)
            return(ret)
        }
        ### conditional models
        mcov <- coef(object, newdata = newdata, type = "Sigma")
        msd <- sqrt(mvtnorm::diagonals(mcov)[margins,])
        if (length(unique(msd)) == 1L && 
            !"bscaling" %in% names(tmp$model$model)) { ### no stram model
            cf <- cf / unique(msd)
            coef(tmp) <- cf
            ret <- predict(tmp, newdata = newdata, type = type, log = log, ...)
            return(ret)
        }
        type <- match.arg(type)
        tr <- predict(tmp, newdata = newdata, type = "trafo", ...) 
        msd <- matrix(msd, nrow = nrow(tr), ncol = ncol(tr), byrow = TRUE)
        tr <- tr / msd
        switch(type, "trafo" = return(tr),
                     "distribution" = return(pnorm(tr, log.p = log)),
                     "survivor" = return(pnorm(tr, log.p = log, 
                                               lower.tail = FALSE)),
                     "density" = {
                         dx <- 1
                         names(dx) <- variable.names(tmp)[1L]
                         dtr <- predict(tmp, newdata = newdata, 
                                        type = "trafo", deriv = dx, ...)
                         ret <- dnorm(tr, log = TRUE) - .log(msd) + .log(dtr)
                         if (log) return(ret)
                         return(exp(ret))
                     },
                     stop("not yet implemented"))
    }

    type <- match.arg(type)
    ### don't feed ...
    z <- .mget(object$models, margins, parm = coef(object, type = "all"),
               newdata = newdata, what = "z")
    z <- .rbind(z)

    if (type == "trafo") {
        stopifnot(!log)
        L <- coef(object, newdata = newdata, type = "Lambda")
        if (length(margins) != J) 
            L <- marg_mvnorm(invchol = L, which = margins)$invchol
        return(Mult(L, z))
    }
    if (type == "distribution") {
        lower <- matrix(-Inf, ncol = ncol(z), nrow = nrow(z))
        upper <- z
        Linv <- coef(object, newdata = newdata, type = "Lambdainv")
        if (length(margins) != J) 
            Linv <- marg_mvnorm(chol = Linv, which = margins)$chol
        a <- args
        a$lower <- lower
        a$upper <- upper
        a$logLik <- FALSE
        a$chol <- Linv
        ret <- do.call("lpmvnorm", a)
        if (log) return(ret)
        return(exp(ret))
    }
    if (type == "survivor") {
        lower <- z 
        upper <- matrix(Inf, ncol = ncol(z), nrow = nrow(z))
        Linv <- coef(object, newdata = newdata, type = "Lambdainv")
        if (length(margins) != J) 
            Linv <- marg_mvnorm(chol = Linv, which = margins)$chol
        a <- args
        a$lower <- lower
        a$upper <- upper
        a$logLik <- FALSE
        a$chol <- Linv
        ret <- do.call("lpmvnorm", a)
        if (log) return(ret)
        return(exp(ret))
    }
    stopifnot(type == "density")
    stopifnot(all(object$models$cresp))
    zprime <- .mget(object$models, margins, parm = coef(object, type = "all"),
                    newdata = newdata, what = "zprime")
    if (length(margins) > 1L) {
        zprime <- .rbind(zprime)
    } else {
        zprime <- matrix(zprime, nrow = 1)
    }
    L <- coef(object, newdata = newdata, type = "Lambda")
        if (length(margins) != J) 
            L <- marg_mvnorm(invchol = L, which = margins)$invchol
    ret <- ldmvnorm(obs = z, invchol = L, logLik = FALSE)
    ret <- ret + colSums(.log(zprime))
    if (log) return(ret)
    return(exp(ret))
}

simulate.mmlt <- function(object, nsim = 1L, seed = NULL, newdata, K = 50, 
                          ...) {

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
            tmp <- object$models$models[[j]]
            q <- mkgrid(tmp, n = K)[[1L]]
            cf <- coef(tmp, fixed = TRUE)
            ncf <- names(cf)
            names(cf) <- paste(variable.names(tmp)[1L], names(cf), sep = ".")
            cfm <- object$models$parm(coef(object, fixed = TRUE))[[j]]
            cf[names(cfm)] <- cfm
            names(cf) <- ncf
            coef(tmp) <- cf
            pr <- predict(tmp, newdata = newdata, type = "trafo", q = q)
            if (!is.matrix(pr)) 
                pr <- matrix(pr, nrow = length(pr), ncol = NROW(newdata))
            ret[,j] <- as.double(mlt:::.invf(tmp, f = t(pr), q = q, 
                                             z = t(Ztilde[j,,drop = FALSE])))
        }
    } else {
        Ztilde <- pnorm(Ztilde, log.p = TRUE)
        for (j in 1:J) {
            tmp <- object$models$models[[j]]
            q <- mkgrid(tmp, n = K)[[1L]]
            cf <- coef(tmp, fixed = TRUE)
            ncf <- names(cf)
            names(cf) <- paste(variable.names(tmp)[1L], names(cf), sep = ".")
            cfm <- object$models$parm(coef(object, fixed = TRUE))[[j]]
            cf[names(cfm)] <- cfm
            names(cf) <- ncf
            coef(tmp) <- cf
            pr <- predict(tmp, newdata = newdata, type = "logdistribution", q = q)
            if (!is.matrix(pr)) 
                pr <- matrix(pr, nrow = length(pr), ncol = NROW(newdata))
            ret[,j] <- as.double(mlt:::.invf(tmp, f = t(pr), q = q, 
                                             z = t(Ztilde[j,,drop = FALSE])))
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

    Linv <- coef(object, newdata = newdata, type = "Lambdainv")
    Linv <- as.array(Linv)[,,1]
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
        tmp <- object$models$models[[j]]
        prb <- tmp$todistr$p(a[,j])
        cf <- coef(tmp)
        cf[] <- object$models$parm(coef(object))[[j]]
        coef(tmp) <- cf
        predict(tmp, newdata = nd, type = "quantile", prob = prb)
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

    ret <- do.call("expand.grid", lapply(object$models$models, 
                                         function(x) mkgrid(x, n = K)[[1L]]))
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
