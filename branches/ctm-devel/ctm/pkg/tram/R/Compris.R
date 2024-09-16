
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

    x <- unclass(object)
    ev <- attr(object, "states")
    J <- length(ev)
    st <- x[, "status"]
    st <- factor(st, levels = 0:J, labels = c("rc", ev))

    if (attr(object, "type") == "mright") {
        ret <- lapply(1:J, function(j) {
            tm <- x[, "time"]
            return(Surv(time = tm, event = st == ev[j]))
        })
        names(ret) <- paste0("Event_", ev)
        return(ret)
    }

    stopifnot(attr(object, "type") == "mcounting")


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
                    NPlogLik = FALSE, theta = NULL,
                    optim = mltoptim(auglag = list(maxtry = 5)), args = list(seed = 1, M = 1000), 
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

    stopifnot(!NPlogLik)

    tmp <- mf
    y <- .Surv2Survs(y)
    tmp[names(y)] <- y

    m <- lapply(1:length(y), function(j) {
        fmj <- formula
        fmj[[2]] <- as.name(names(y)[j])
        if (j == 1)
            return(FUNe(formula = fmj, data = tmp, ...))
        return(FUNc(formula = fmj, data = tmp))
    })
    names(m) <- names(y)
    m$data <- tmp
    m$optim <- optim
    m$theta <- theta
    m$scale <- scale
    m$args <- args
    ret <- do.call("Mmlt", m)
    ret$call <- call
    class(ret) <- c("Compris", class(ret))
    return(ret)
}
