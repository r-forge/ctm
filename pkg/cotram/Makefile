
PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)

PKGFiles := $(shell ls)

update: install .Rout.save example.Rout.save

update_dep: 
	echo "install.packages(c('basefun', 'variables', 'mlt', 'tram'), repos = 'https://cran.r-project.org')" | R --vanilla

compactPDF:
	cd ..;\
	compactPDF $(PKGNAME)_$(PKGVERS).tar.gz

../$(PKGNAME)_$(PKGVERS).tar.gz:	$(PKGFiles) .Rbuildignore
	cd ..; R CMD build --compact-vignettes=both $(PKGSRC)
	# make compactPDF
	
build:	../$(PKGNAME)_$(PKGVERS).tar.gz

install:	../$(PKGNAME)_$(PKGVERS).tar.gz
	cd ..;\
	R CMD INSTALL $(PKGNAME)_$(PKGVERS).tar.gz

check: install
	cd ..;\
	R CMD check --as-cran $(PKGNAME)_$(PKGVERS).tar.gz
	
cotram.Rcheck: check
	
checkdevel: build
	cd ..; Rdev CMD check --as-cran $(PKGNAME)_$(PKGVERS).tar.gz --as-cran
	## https://www.r-bloggers.com/2015/07/r-devel-in-parallel-to-regular-r-installation/

testfiles = confband-Ex.R # mcotram-Ex.R
	
.Rout.save:	$(PKGFiles) .Rbuildignore cotram.Rcheck 
	cp ../$(PKGNAME).Rcheck/$(PKGNAME)-Ex.Rout tests/Examples/$(PKGNAME)-Ex.Rout.save
	for file in $(testfiles); do \
		cp ../$(PKGNAME).Rcheck/tests/$${file}out tests/$${file%.R}.Rout.save; \
	done

clean:
	cd ..;\
	rm -rf ../$(PKGNAME).Rcheck/
	rm -rf ../$(PKGNAME)_$(PKGVERS).tar.gz
	./cleanup

