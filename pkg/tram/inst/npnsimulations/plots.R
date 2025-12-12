
library("lattice")

load("summary.rda")

lev <- levels(args$method)
nlev <- expression(ML, MPL, cML, "\u2113"[2](theta, lambda), tilde("\u2113")[2](theta(vartheta), lambda), tilde("\u2113")[2](vartheta, lambda))

args$r <- args$rho <- round(with(args, -lambda / sqrt(1 + lambda^2)), 1)
args$ncat1 <- factor(args$ncat1, levels = c(2, 5, Inf), labels = c("Binary", "Ordinal", "Cont."))
args$ncat2 <- factor(args$ncat2, levels = c(2, 5, Inf), labels = c("Binary", "Ordinal", "Cont."))
args$ncat <- with(args, interaction(ncat1, ncat2, sep = " / "))[,drop = TRUE]
args$n <- args$N
v <- c("N", "rho")#, "ncat")
args[,v] <- lapply(v, function(var) {
    lev <- sort(unique(args[[var]]))
    lab <- paste(var, lev, sep = "=")
    factor(args[[var]], levels = lev, labels = lab)
})

failed <- is.na(args$Est)
xtabs(~ method + N + rho + ncat, data = args[failed,,drop = FALSE])

cairo_pdf("figures.pdf", width = 10, height = 10)

for (r in levels(args$rho)) {
  print(r)
  tmp <- args[args$rho == r,]# & args$method != "mvord",]
  plot(bwplot(Est ~ method | N + ncat, data = tmp, ylim = c(-1, 1),
              varwidth = TRUE,
       panel = function(..., subscripts = subscripts) {
           r <- unique(tmp[subscripts, "r"])
           panel.abline(h = r, col = "darkgrey", lty = 3)
           panel.bwplot(...)
       }, main = parse(text = gsub("=", "==", r)), ylab = expression(hat(rho)),
       scales = list(y = list(rot = 0, at = c(-1, 0, 1), labels = c(-1, 0, 1)), 
                     x = list(at = 1:6, labels = nlev, rot = 45))))
  plot(bwplot(abs(SE) ~ method | N + ncat, data = tmp, ylim = c(0, .5),
            subset = abs(SE) < 100, varwidth = TRUE,
            panel = function(..., subscripts = subscripts) {
               n <- unique(tmp[subscripts, "n"])
               r <- unique(tmp[subscripts, "r"])
               est <- tmp[subscripts, "Est"]
               m <- tmp[subscripts, "method"]
               sdm <- tapply(est, m, sd, na.rm = TRUE)
               panel.abline(h = sqrt((1 - r^2)^2 / n), col = "darkgrey", lty = 3)
               panel.bwplot(...)
               panel.points(x = 1:length(sdm), y = sdm, col = "red", 
                            pch = 17)
       }, main = parse(text = gsub("=", "==", r)), ylab = expression(SE(hat(rho))),
       scales = list(y = list(rot = 0), #relation = "free", log = FALSE),
                     x = list(at = 1:6, labels = nlev, rot = 45))))
  plot(bwplot(Time ~ method | N + ncat, data = tmp, varwidth = TRUE,
        main = parse(text = gsub("=", "==", r)), scales = list(y = list(log = TRUE),
                               x = list(at = 1:6, labels = nlev, rot = 45))
       # , scales = list(y = list(relation = "free")
       ))
}

dev.off()
