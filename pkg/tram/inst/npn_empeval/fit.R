
library("tram")
library("mvtnorm")
library("qrng")

fit <- function(data, mfun = BoxCox, as.R.interval = FALSE, as.R.ordered = FALSE,
                optim = mltoptim(trace = FALSE), 
                M = 500,
                conditional = FALSE, 
                domargins = TRUE, sequentialfit = FALSE, se = FALSE, seed = NULL, ...) { 

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

    marg <- lapply(1:ncol(data), function(j) {
        if (is.factor(data[[j]])) {
            fm <- as.formula(paste(colnames(data)[j], " ~ 1"))
            return(Polr(fm, data = data, method = "probit", ...))
        }
        fm <- as.formula(paste(colnames(data)[j], " ~ 1"))
        if (as.R.interval)
            fm <- as.formula(paste("R(", colnames(data)[j], ", 
                                   as.R.interval = TRUE) ~ 1"))
        if (as.R.ordered)
            fm <- as.formula(paste("R(", colnames(data)[j], ", 
                                   as.R.ordered = TRUE) ~ 1"))
        return(mfun(fm, data = data, ...))
    })

    m <- marg
    m$formula <- ~ 1
    m$data <- data
    m$conditional <- conditional
    m$domargins <- domargins
    m$sequentialfit <- sequentialfit
    m$args <- list(seed = 1, type = c("ghalton"), M = M)
    ### might want to switch to nloptr if hessian is not necessary
    # if (!se) optim <- optim["nloptr"]
    m$optim <- optim
    ret <- do.call("Mmlt", m)
    L <- coef(ret, type = "Lambdapar")
    if (!se) return(L)

    mp <- sum(sapply(lapply(marg, function(x) coef(as.mlt(x), fixed = FALSE)), length))
    se <- sqrt(diag(vcov(ret)))
    se <- se[-(1:mp)]
    se <- ltMatrices(se, byrow = TRUE)
    return(list(L = L, se = se))
}
