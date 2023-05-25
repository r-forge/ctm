
PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)

update: install .Rout.save example.Rout.save

compactPDF:
	cd ..;\
	compactPDF $(PKGNAME)_$(PKGVERS).tar.gz

../$(PKGNAME)_$(PKGVERS).tar.gz:
	cd ..; R CMD build $(PKGSRC)
	make compactPDF
	
build:	../$(PKGNAME)_$(PKGVERS).tar.gz
	make ../$(PKGNAME)_$(PKGVERS).tar.gz

install: build
	cd ..;\
	R CMD INSTALL $(PKGNAME)_$(PKGVERS).tar.gz

check: build
	cd ..;\
	R CMD check $(PKGNAME)_$(PKGVERS).tar.gz --as-cran
	
checkdevel: build
	cd ..; Rdev CMD check --as-cran $(PKGNAME)_$(PKGVERS).tar.gz --as-cran
	## https://www.r-bloggers.com/2015/07/r-devel-in-parallel-to-regular-r-installation/

testfiles = confband-Ex.R mcotram-Ex.R

installall:
	make -C ../ installall

.Rout.save:
	for file in $(testfiles); do \
		cd tests; R CMD BATCH --vanilla $${file} $${file%.R}.Rout.save; \
	done
	
example.Rout.save:
	cd ..; R CMD check $(PKGNAME)_$(PKGVERS).tar.gz
	cd ..; cp $(PKGNAME).Rcheck/$(PKGNAME)-Ex.Rout $(PKGNAME)/tests/Examples/$(PKGNAME)-Ex.Rout.save

clean:
	cd ..;\
	# $(RM) -r $(PKGNAME).Rcheck/
	rm -rf ../$(PKGNAME).Rcheck/
	rm -rf ../$(PKGNAME)_$(PKGVERS).tar.gz

