
source("setup.R")
library("copula")
nc <- normalCopula(dim = 2, dispstr = "un")

set.seed(2908)

load("data.rda")

ret <- matrix(NA, nrow = nrow(args), ncol = 3)
colnames(ret) <- c("Est", "SE", "Time")

run <- function(i) {
    print(args[i,])

    ret <- matrix(NA, nrow = 1, ncol = 3)
    colnames(ret) <- c("Est", "SE", "Time")

    ### only continuous data
    if (any(sapply(d[[i]], is.factor))) return(ret)
    
    ret[1, "Time"] <- system.time(mm <- try(fitCopula(nc, data = pobs(d[[i]]), method = "ml")))["user.self"]
    if (inherits(mm, "try-error")) {
        ret[1, "Time"] <- NA
        return(ret)
    }

    ret[1,"Est"] <- coef(mm)
    ret[1,"SE"] <- sqrt(vcov(mm))
    return(ret)
}

ret <- do.call("rbind", mclapply(1:nrow(args), run, mc.cores = MC))

save(ret, file = "ret_eff_copula_ml.rda")
