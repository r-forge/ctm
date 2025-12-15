
load("data.rda")

args$Est <- args$SE <- args$Time <- 0
vars <- c("Est", "SE", "Time")

method <- factor(c("P_P", "NP_P", "NP_NP", "copula", "copula_ml", "mvord"))

tmp <- expand.grid(i = 1:nrow(args), method = method)

args <- args[tmp$i,]
args$method <- tmp$method

load("ret_eff_P_P.rda")
args[args$method == "P_P", vars] <- as.numeric(ret)

load("ret_eff_NP_P.rda")
args[args$method == "NP_P", vars] <- as.numeric(ret)

load("ret_eff_NP_NP.rda")
args[args$method == "NP_NP", vars] <- as.numeric(ret)

load("ret_eff_copula.rda")
args[args$method == "copula", vars] <- as.numeric(ret)

load("ret_eff_copula_ml.rda")
args[args$method == "copula_ml", vars] <- as.numeric(ret)

load("ret_eff_mvord.rda")
args[args$method == "mvord", vars] <- as.numeric(ret)

save(args, file = "summary.rda")
