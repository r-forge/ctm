
data.rda:	mkdata.R
	R CMD BATCH mkdata.R

R = $(shell ls  eff*.R)
Rout = $(R:.R=.Rout)

.SUFFIXES:
.SUFFIXES: .R .Rout

.R.Rout:	
	R CMD BATCH $< 

summary.rda:	summary.R $(Rout)
	R CMD BATCH summary.R

figures.pdf:	summary.rda plots.R
	R CMD BATCH plots.R

all: data.rda $(Rout) summary.rda figures.pdf
	
clean:	
	rm -rf *.rda
	rm -rf *.Rout
	rm -rf *~
	rm -rf figures.pdf

