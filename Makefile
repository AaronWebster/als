ALS_VERSION := 5.0.0

MKOCTFILE := CXXFLAGS=-g mkoctfile -v

OCTAVE := octave -qf

NMPC_DISTFILE := mpc-tools-$(ALS_VERSION)

INSTALL = /usr/bin/install -c
INSTALL_PROGRAM = $(INSTALL)
INSTALL_SCRIPT = $(INSTALL})
INSTALL_DATA = $(INSTALL) -m 644

OCT_INSTDIR := $(DESTIDR)$(shell octave-config -p LOCALAPIOCTFILEDIR)/als

M_INSTDIR := $(DESTDIR)$(shell octave-config -p LOCALAPIFCNFILEDIR)/als

M_SRC := \
	als_diag.m \
	golden_section_Q_mrQ.m \
	obj_tot.m \
	obj_ls.m \
	symtran.m \
	symtranT.m \
	als_sdp_mrQ.m \
	minrank_ex.m \
	sdp_QR_mrQ.m \
	celldiag.m \
	simulate_data8_diag.m \
	ECM_iter.m \
	motivating_example.m \
	simulate_data8.m \
	ltv_als.m \
	case1_als.m \
	cvode_sens24.m \
	als_condno.m \
	autocov_calc.m \
	autocov_calc_thry.m \
	comm_mat.m \
	xdata_case1.mat \
	ydata_case1.mat

CXX_SRC := 

CXX_OBJ := $(CXX_SRC:.cc=.o)
CXX_OCT := $(CXX_SRC:.cc=.oct)

EXTRAS := 

DISTFILES := $(M_SRC) $(CXX_SRC) $(EXTRAS)

DISTSUBDIRS := test

all: 
.PHONY: all

dist:
	$(eval TMP := $(shell mktemp -d))
	@mkdir $(TMP)/als
	@mkdir $(TMP)/als/src
	@mkdir $(TMP)/als/inst
	echo "Name: als" > $(TMP)/als/DESCRIPTION
	echo "Version: 5.0.0" >> $(TMP)/als/DESCRIPTION
	echo "Date: 2014-05-12" >> $(TMP)/als/DESCRIPTION
	echo "Author: Fernando V. Lima and Murali R. Rajamani" >> $(TMP)/als/DESCRIPTION
	echo "Maintainer: Fernando V. Lima (flima@bevo.che.wisc.edu)" >> $(TMP)/als/DESCRIPTION
	echo "Title: ALS Toolbox" >> $(TMP)/als/DESCRIPTION
	echo "Description: test" >> $(TMP)/als/DESCRIPTION
	echo "Autoload: yes" >> $(TMP)/als/DESCRIPTION
	echo "Categories: system identification tools" >> $(TMP)/als/DESCRIPTION
	echo "SystemRequirements: none" >> $(TMP)/als/DESCRIPTION
	cp $(M_SRC) $(TMP)/als/inst/
	cp Makefile $(TMP)/als/src/
	cp README.md $(TMP)/als/README
	cp LICENSE $(TMP)/als/COPYING
	cd $(TMP) && tar -czvf als.tar.gz als/
	mv $(TMP)/als.tar.gz .
	rm -rf $(TMP)
.PHONY: dist

clean:
	rm -f als.tar.gz
.PHONY: clean
