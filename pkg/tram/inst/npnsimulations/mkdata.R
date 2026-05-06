
set.seed(29)
source("dgp.R")

Nsim <- 100

r <- 0:9 / 10
lambda <- -1 * sign(r) * sqrt(r^2 / (1 - r^2))

args <- expand.grid(N = c(10, 30, 50), lambda = lambda, ncat1 = c(Inf, 5, 2), ncat2 = c(Inf, 5, 2), Nsim = 1:Nsim)
args <- subset(args, ncat1 >= ncat2)
args$sim <- 1:nrow(args)
args$seed <- floor(runif(nrow(args)) * 10^6)
args$id <- 1:nrow(args)
d <- vector(mode = "list", length = nrow(args))

for (i in 1:nrow(args))
    d[[i]] <- dgp(n = args[i, "N"], df = 2, lambda = args[i, "lambda"], 
                  ncat = c(args[i, "ncat1"], args[i, "ncat2"]))

print(length(d))

save(args, d, file = "data.rda")
