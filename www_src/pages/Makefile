
Rmd = $(shell ls *.Rmd)
md = $(Rmd:.Rmd=.md)

.SUFFIXES:
.SUFFIXES: .Rmd .md

.Rmd.md:
	echo "download.file(\"https://user.math.uzh.ch/hothorn/bib/TH.bib\", dest = \"TH.bib\"); library(\"knitr\"); knit(\"$<\");" | R --vanilla 
	
md: $(md)

clean:
	rm -rf $(md)
	rm -rf TH.bib
	rm -rf *~
	
