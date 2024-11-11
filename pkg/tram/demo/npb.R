
### Reproducibility material for
###
### Nonparanormal Modeling Framework for
### Prognostic Biomarker Assessment
###
### ALS study
pdf("npb.pdf")
source("https://gitlab.com/asewak/npb/-/raw/main/als.R?ref_type=heads&inline=false", echo = TRUE)
dev.off()
###
### Simulation; takes a while
if (FALSE) 
    source("https://gitlab.com/asewak/npb/-/raw/main/sim.R?ref_type=heads&inline=false", echo = TRUE)

sessionInfo()
