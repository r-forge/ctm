
.mrightORmcounting <- function(object) {
    if (!inherits(object, "Surv")) return(FALSE)
    type <- attr(object, "type")
    if (!(type %in% c("mright", "mcounting")))
        return(FALSE)
    st <- unclass(object)[, "status"]
    if (all(unique(st) %in% c(0, 1))) return(FALSE)
    return(TRUE)
}

.mright2NP <- function(object) {
    stopifnot(inherits(object, "Surv"))
    stopifnot(attr(object, "type") == "mright")
    stopifnot(!is.null(ev <- attr(object, "states")))

    time <- object[,"time"]
    status <- object[, "status"]
    idx <- status == 0
    TE <- time[idx]
    J <- 1:length(ev)
    ret <- lapply(J, function(j) {
        ret <- as.Surv(R(Surv(time, event = status == j), 
                         as.R.interval = TRUE))
        ret[idx,"time1"] <- TE
        ret
    })
    for (j in J) {
        idx <- status == j
        TE <- time[idx]
        for (k in J[J != j])
            ret[[k]][idx,"time1"] <- TE
    }
    names(ret) <- paste0("Event_", ev)
    return(ret)
}

.Surv2Survs <- function(object) {

    stopifnot(attr(object, "type") == "mcounting")

    x <- unclass(object)
    ev <- attr(object, "states")
    J <- length(ev)
    st <- x[, "status"]
    st <- factor(st, levels = 0:J, labels = c("rc", ev))

    NAstop <- is.na(x[, "stop"])

    if (!any(NAstop[st != "rc"])) {
        ret <- lapply(1:J, function(j) {
            tm <- x[, "start"]
            tm2 <- x[, "stop"]
            ### rc
            tm2[st == "rc"] <- Inf
            ### other event
            OE <- st %in% ev[-j]
            tm[OE] <- tm2[OE]
            tm2[OE] <- Inf
            return(Surv(time = tm, time2 = tm2, type = "interval2"))
        })
        names(ret) <- paste0("Event_", ev)
        return(ret)
    }
}

Compris <- function(formula, data, subset, weights, na.action, offset,
                    primary = c("Coxph", "Colr", "BoxCox"),
                    competing = switch(primary, "Coxph" = "weibull", 
                                                "Colr" = "loglogistic", 
                                                "BoxCox" = "lognormal"),
                    NPlogLik = FALSE, 
                    optim = mmltoptim(), args = list(seed = 1, M = 1000), 
                    scale = FALSE, tol = .001, ...)
{

    call <- match.call()

    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action", "weights", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(model.frame)
    mf <- eval(mf, parent.frame())

    y <- model.response(mf)
    w <- as.vector(model.weights(mf))
    if (!is.null(w) && !is.numeric(w)) 
        stop("'weights' must be a numeric vector")
    offset <- model.offset(mf)

    primary <- match.arg(primary)
    FUNe <- function(...) do.call(primary, list(...))
    FUNc <- function(...) Survreg(..., dist = competing)

    stopifnot(inherits(y, "Surv"))
    tm <- y[,1]
    st <- attr(y, "states")
    ev <- factor(y[,ncol(y)], levels = 0:length(st), labels = c("admin_cens", st))
    J <- nlevels(ev) - 1L
    if (nlevels(ev) == 2)
        return(FUNe(formula = formula, data = mf, order = order[J], ...))

    if (!(attr(y, "type") == "mright" && J == 2 && !NPlogLik)) {
        ### all observations are intervals, maximise J-dim probabilities
        if (attr(y, "type") == "mcounting") {
            tmp <- mf
            y <- .Surv2Survs(y)
            tmp[names(y)] <- y
        } else {
            ### J > 2 && attr(y, "type") == "mright" || NPlogLik
            tmp <- mf
            y <- .mright2NP(y)
            tmp[names(y)] <- y
        }

        m <- lapply(1:length(y), function(j) {
            fmj <- formula
            fmj[[2]] <- as.name(names(y)[j])
            if (j == 1)
                return(FUNe(formula = fmj, data = tmp, ...))
            return(FUNc(formula = fmj, data = tmp))
        })
        names(m) <- names(y)

        m$optim <- optim
        m$scale <- scale
        m$args <- args
        ret <- do.call("mmlt", m)
        ret$call <- call
        class(ret) <- c("Compris", class(ret))
        return(ret)
    }

    ### mright and J == 2 with tricks

    stopifnot(attr(y, "type") == "mright")
    if (scale)
        warning("scaling of coefficients not yet implemened")

    tmp <- mf
    tmp[[st[1]]] <- with(tmp, Surv(y[,1], event = ev == st[1]))
    fme <- formula
    fme[[2]] <- as.name(st[1])
    me <- as.mlt(FUNe(fme, data = tmp, ...))
    cf_me <- coef(me)
    tmp[[st[2]]] <- with(tmp, R(Surv(y[,1], event = ev == st[2])))
    fmc <- formula
    fmc[[2]] <- as.name(st[2]) ### Surv makes assumptions
                               ### about bounds and support
    mc <- as.mlt(FUNc(fmc, data = tmp))
    cf_mc <- coef(mc)

    ### construct log-likelihood as sum of three terms:
    # A) administrative right-censoring
    rc <- which(y[,2] == 0)
    if (length(rc) > 0) {
        d_rc <- mf[rc,,drop = FALSE]
        d_rc[[st[1]]] <- d_rc[[st[2]]] <- Surv(d_rc[[1]][,1], event = rep(0, length(rc)))
    }
    # B) event of interest
    ei <- which(y[,2] == 1)
    d_ei <- mf[ei,,drop = FALSE]
    d_ei[[st[1]]] <- Surv(d_ei[[1]][,1], event = rep(1, length(ei)))
    d_ei[[st[2]]] <- Surv(d_ei[[1]][,1], event = rep(0, length(ei)))
    # C) competing event
    ce <- which(y[,2] == 2)
    d_ce <- mf[ce,,drop = FALSE]
    d_ce[[st[1]]] <- Surv(d_ce[[1]][,1], event = rep(0, length(ce)))
    d_ce[[st[2]]] <- Surv(d_ce[[1]][,1], event = rep(1, length(ce)))

    ### set-up models for each of these subsets
    SW <- suppressWarnings
    if (length(rc) > 0)
        SW(me_rc <- mlt(me$model, data = d_rc, theta = cf_me, dofit = FALSE))
    SW(me_ei <- mlt(me$model, data = d_ei, theta = cf_me, dofit = FALSE))
    SW(me_ce <- mlt(me$model, data = d_ce, theta = cf_me, dofit = FALSE))
    if (length(rc) > 0)
        SW(mc_rc <- mlt(mc$model, data = d_rc, theta = cf_mc, dofit = FALSE))
    SW(mc_ei <- mlt(mc$model, data = d_ei, theta = cf_mc, dofit = FALSE))
    SW(mc_ce <- mlt(mc$model, data = d_ce, theta = cf_mc, dofit = FALSE))

    ### set-up likelihood terms
    # A)
    if (length(rc) > 0) {
        mm_rc <- mmlt(me_rc, mc_rc, dofit = FALSE, args = args)
    }
    # B)
    mm_ei <- mmlt(me_ei, mc_ei, dofit = FALSE)
    # C) NOTE that the model order is different here!
    mm_ce <- mmlt(mc_ce, me_ce, dofit = FALSE)

    ro <- function(parm, reverse = FALSE) {
        if (!reverse) {
            p_ei <- parm[1:length(cf_me)]
            p_cp <- parm[length(cf_me) + 1:length(cf_mc)]
            p_cr <- parm[-(1:(length(cf_me) + length(cf_mc)))]
            return(c(p_cp, p_ei, p_cr))
        }
        p_cp <- parm[1:length(cf_mc)]
        p_ei <- parm[length(cf_mc) + 1:length(cf_me)]
        p_cr <- parm[-(1:(length(cf_me) + length(cf_mc)))]
        return(c(p_ei, p_cp, p_cr))
    }

    ll <- function(parm) {
        ret <- 0
        if (length(rc) > 0) ret <- ret + mm_rc$ll(parm)
        ret <- ret + mm_ei$ll(parm)
        ret <- ret + mm_ce$ll(ro(parm))
        ret
    }

    sc <- function(parm) {
        ret <- 0
        if (length(rc) > 0) ret <- ret + mm_rc$sc(parm)
        ret <- ret + mm_ei$sc(parm)
        ret <- ret + ro(mm_ce$sc(ro(parm)), reverse = TRUE)
        ret
    }

    theta <- c(coef(me), coef(mc), lambda = 0)

    for (i in 1:length(optim)) {
        op <- try(optim[[i]](theta, ll, sc, mm_ei$ui, mm_ei$ci))
        ## <FIXME> maybe add more informative error message </FIXME>
        if (inherits(op, "try-error")) op <- list(convergence = 1)
        if (op$convergence == 0) break()
    }
    names(op$par) <- names(theta)    

    ### fake global mmlt model: this just takes ret$par
    ### and doesn't optimise further
    ret <- mmlt(me_ei, mc_ei, data = d_ei, theta = op$par, 
                dofit = FALSE)
    ret$optim_hessian <- NULL
    ret$logLik <- -op$value
    ret$ll <- function(parm, newdata = NULL) {
        if (!is.null(newdata)) warning("Argument newdata ignored")
        return(ll(parm))
    }
    ret$sc <- function(parm, newdata = NULL) {
        if (!is.null(newdata)) warning("Argument newdata ignored")
        return(sc(parm))
    }
    if (!is.null(op$optim_hessian))
        ret$optim_hessian <- op$optim_hessian
    ret$call <- call
    ret$mmlt <- "Competing Risk Regression"
    ret$args <- args
    class(ret) <- c("Compris", class(ret))
    ret
}
