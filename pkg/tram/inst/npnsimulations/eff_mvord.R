
source("setup.R")
library("mvord")

load("data.rda")

run <- function(i) {
    print(args[i,])

    ret <- matrix(NA, nrow = 1, ncol = 3)
    colnames(ret) <- c("Est", "SE", "Time")

    dJ <- sum(sapply(d[[i]], is.factor))
    if (dJ < 2) return(ret)

    tmp <- d[[i]]
    ret[1, "Time"] <- system.time(mm <- try(mvord(formula = MMO2(Y1, Y2) ~ 1, data = tmp, 
                                  control = mvord.control(solver="nlminb"))))["user.self"]
    if (inherits(mm, "try-error")) {
        ret[1, "Time"] <- NA
        return(ret)
    }

    sink(tempfile())
    tmp <- try(summary(mm)$error.structure[1, "Estimate"])
    if (!inherits(tmp, "try-error"))
        ret[1,"Est"] <- tmp 
    tmp <- try(summary(mm)$error.structure[1,"Std. Error"])
    if (!inherits(tmp, "try-error"))
        ret[1,"SE"] <- tmp
    sink()
    return(ret)
}

ret <- do.call("rbind", mclapply(1:nrow(args), run, mc.cores = MC))

save(ret, file = "ret_eff_mvord.rda")
