# Makefile for R Journal paper 

RNW = $(wildcard *.Rnw)

pres:
	Rscript -e "library(knitr); knit2pdf('$(RNW)')"

clean:
	rm -f *~ *.log *.aux *.out Rplots.pdf \
	  *.synctex.gz .RData .Rhistory *-concordance*.tex \
	  *.bbl *.nav *.toc *.snm *.vrb \
	  *.blg *.brf $(RNW:.Rnw=.tex)

cleanall: clean
	rm -rf cache/ files/ figure/

all: pres clean

