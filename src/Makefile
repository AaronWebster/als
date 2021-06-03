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

PKG_SRC := COPYING \
	README \
	DESCRIPTION

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
	xdata_case1.dat \
	ydata_case1.dat

CXX_SRC := 

CXX_OBJ := $(CXX_SRC:.cc=.o)
CXX_OCT := $(CXX_SRC:.cc=.oct)

EXTRAS := 

DISTFILES := $(M_SRC) $(CXX_SRC) $(EXTRAS)

DISTSUBDIRS := test

all: 
.PHONY: all

dist:
	if test -d als; then echo; else mkdir als; fi
	if test -d als/src; then echo; else mkdir als/src; fi
	if test -d als/inst; then echo; else mkdir als/inst; fi
	cp $(M_SRC) als/inst/
	cp Makefile als/src/
	cp $(PKG_SRC) als/
	tar -czvf als.tar.gz als/
	rm -rf als
.PHONY: dist

clean:
	rm -f als.tar.gz
	rm -f *.*~
	rm -f *~
.PHONY: clean
