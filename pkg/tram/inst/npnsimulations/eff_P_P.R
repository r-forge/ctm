
source("setup.R")
library("tram")
source("fit.R")

load("data.rda")

run <- function(i) {
    print(args[i,])

    ret <- matrix(NA, nrow = 1, ncol = 3)
    colnames(ret) <- c("Est", "SE", "Time")

    ### at least one continuous variable
    if (all(sapply(d[[i]], is.factor))) return(ret)

    ret[1, "Time"] <- system.time(mm <- try(fit(d[[i]], se = TRUE)))["user.self"]
    if (inherits(mm, "try-error")) {
        ret[1, "Time"] <- NA
        return(ret)
    }

    ret[1,"Est"] <- as.array(invchol2cor(mm$L))[1,2,1]
    lam <- unclass(mm$L)
    lams2 <- unclass(mm$se)
    ret[1,"SE"] <- lams2 * abs(lam^2 / ((lam^2 + 1)^(3/2)) - 1 / sqrt(lam^2 + 1))
    return(ret)
}

ret <- do.call("rbind", mclapply(1:nrow(args), run, mc.cores = MC))

save(ret, file = "ret_eff_P_P.rda")
