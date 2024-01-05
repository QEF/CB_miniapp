# Copyright (C) 2001-2016 Quantum ESPRESSO group
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.

include make.inc

default :
	@echo 'to install this , type at the shell prompt:'
	@echo '  ./configure [--prefix=]'
	@echo '  make [-j] target'

###########################################################
# Main targets
###########################################################

# The syntax "( cd PW ; $(MAKE) TLDEPS= all || exit 1)" below
# guarantees that error code 1 is returned in case of error and make stops
# If "|| exit 1" is not present, the error code from make in subdirectories
# is not returned and make goes on even if compilation has failed

cb : bindir libfft libdavid libdavid_rci libcg libppcg libparo libdense libla liblapack libutil
	if test -d CB_toy_code ; then \
	( cd CB_toy_code ; $(MAKE) TLDEPS= all || exit 1) ; fi

###########################################################
# Auxiliary targets used by main targets:
# compile modules, libraries, directory for binaries, etc
###########################################################

libdavid_rci : touch-dummy libla clib libutil
	( cd KS_Solvers/Davidson_RCI ; $(MAKE) TLDEPS= all || exit 1 )

libdavid : touch-dummy libla clib libutil
	( cd KS_Solvers/Davidson ; $(MAKE) TLDEPS= all || exit 1 )

libcg : touch-dummy libla clib libutil
	( cd KS_Solvers/CG ; $(MAKE) TLDEPS= all || exit 1 )

libppcg : touch-dummy libla clib libutil
	( cd KS_Solvers/PPCG ; $(MAKE) TLDEPS= all || exit 1 )

libparo : touch-dummy libla clib libutil
	( cd KS_Solvers/ParO ; $(MAKE) TLDEPS= all || exit 1 )

libdense : touch-dummy libla clib libutil
	( cd KS_Solvers/DENSE ; $(MAKE) TLDEPS= all || exit 1 )

libla : touch-dummy liblapack libutil
	( cd LAXlib ; $(MAKE) TLDEPS= all || exit 1 )

libfft : touch-dummy
	( cd FFTXlib ; $(MAKE) TLDEPS= all || exit 1 )

libutil : touch-dummy
	( cd UtilXlib ; $(MAKE) TLDEPS= all || exit 1 )

clib : touch-dummy
	( cd clib ; $(MAKE) TLDEPS= all || exit 1 )

bindir :
	test -d bin || mkdir bin

#############################################################
# Targets for external libraries
############################################################

libblas : touch-dummy
	cd install ; $(MAKE) -f extlibs_makefile $@

liblapack: touch-dummy
	cd install ; $(MAKE) -f extlibs_makefile $@


#########################################################
# plugins
#########################################################

touch-dummy :
	$(dummy-variable)

#########################################################
# "make links" produces links to all executables in bin/
#########################################################

# Contains workaround for name conflicts (dos.x and bands.x) with WANT
links : bindir
	( cd bin/ ; \
	rm -f *.x ; \
	for exe in ../*/*/*.x ../*/bin/* ; do \
	    if test ! -L $$exe ; then ln -fs $$exe . ; fi \
	done ; \
	)

#########################################################
# 'make install' works based on --with-prefix
# - If the final directory does not exists it creates it
#########################################################

install : touch-dummy
	@if test -d bin ; then mkdir -p $(PREFIX)/bin ; \
	for x in `find * ! -path "test-suite/*" -name *.x -type f` ; do \
		cp $$x $(PREFIX)/bin/ ; done ; \
	fi
	@echo 'binaries installed in $(PREFIX)/bin'

#########################################################
# Run test-suite for numerical regression testing
# NB: it is assumed that reference outputs have been
#     already computed once (usualy during release)
#########################################################

test-suite: touch-dummy
	( cd install ; $(MAKE) -f plugins_makefile $@ || exit 1 )

#########################################################
# Other targets: clean up
#########################################################

# remove object files and executables
clean :
	touch make.inc
	for dir in \
	    CB_toy_code LAXlib FFTXlib UtilXlib KS_Solvers/Davidson KS_Solvers/Davidson_RCI KS_Solvers/CG KS_Solvers/PPCG KS_Solvers/ParO clib \
	; do \
	    if test -d $$dir ; then \
		( cd $$dir ; \
		$(MAKE) TLDEPS= clean ) \
	    fi \
	done
	- @(cd install ; $(MAKE) -f plugins_makefile clean)
#	- @(cd install ; $(MAKE) -f extlibs_makefile clean)
	- /bin/rm -rf bin/*.x tmp

# remove files produced by "configure" as well
veryclean : clean
	- @(cd install ; $(MAKE) -f plugins_makefile veryclean)
	- @(cd install ; $(MAKE) -f extlibs_makefile veryclean)
	- cd install ; rm -f config.log configure.msg config.status \
		ChangeLog* intel.pcl */intel.pcl
	- rm -rf include/configure.h
	- cd install ; rm -fr autom4te.cache
	- cd install; ./clean.sh ; cd -
	- cd include; ./clean.sh ; cd -
	- rm -f project.tar.gz
	- rm -rf make.inc

# remove everything not in the original distribution
distclean : veryclean
	( cd install ; $(MAKE) -f plugins_makefile $@ || exit 1 )

tar :
	@if test -f eslw_drivers.tar.gz ; then /bin/rm project.tar.gz ; fi
	# do not include unneeded stuff
	find ./ -type f | grep -v -e /.svn/ -e'/\.' -e'\.o$$' -e'\.mod$$'\
		-e /.git/ -e'\.a$$' -e'\.d$$' -e'\.i$$' -e'_tmp\.f90$$' -e'\.x$$' \
		-e'~$$' -e'\./GUI' -e '\./tempdir' | xargs tar rvf eslw_drivers.tar
	gzip eslw_drivers.tar

#########################################################
# Tools for the developers
#########################################################

# NOTICE about "make doc": in order to build the .html and .txt
# documentation in Doc, "tcl", "tcllib", "xsltproc" are needed;
# in order to build the .pdf files in Doc, "pdflatex" is needed;
# in order to build html files for user guide and developer manual,
# "latex2html" and "convert" (from Image-Magick) are needed.
doc : touch-dummy
	if test -d Doc ; then \
	( cd Doc ; $(MAKE) TLDEPS= all ) ; fi
	for dir in */Doc; do \
	( if test -f $$dir/Makefile ; then \
	( cd $$dir; $(MAKE) TLDEPS= all ) ; fi ) ;  done

doc_clean :
	if test -d Doc ; then \
	( cd Doc ; $(MAKE) TLDEPS= clean ) ; fi
	for dir in */Doc; do \
	( if test -f $$dir/Makefile ; then \
	( cd $$dir; $(MAKE) TLDEPS= clean ) ; fi ) ;  done

depend:
	@echo 'Checking dependencies...'
	- ( if test -x install/makedeps.sh ; then install/makedeps.sh ; fi)
# update file containing version number before looking for dependencies
