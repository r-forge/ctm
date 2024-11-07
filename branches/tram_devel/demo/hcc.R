
### Reproducibility material for
###
### Construction and evaluation of optimal diagnostic 
### tests with application to hepatocellular carcinoma diagnosis
###
### Hepatocellular Carcinoma Diagnosis
pdf("hcc.pdf")
source("https://gitlab.com/asewak/optcomb/-/raw/main/hcc.R?ref_type=heads&inline=false", echo = TRUE)
dev.off()
###
### Simulation; takes a while
if (FALSE) 
    source("https://gitlab.com/asewak/optcomb/-/raw/main/sim.R?ref_type=heads&inline=false", echo = TRUE)
