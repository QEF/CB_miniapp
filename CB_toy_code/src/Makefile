# Makefile for PW

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG)../../FFTXlib $(MOD_FLAG)../../LAXlib $(MOD_FLAG) ../../UtilXlib $(MOD_FLAG)../../KS_Solvers/Davidson $(MOD_FLAG)../../KS_Solvers/Davidson_RCI $(MOD_FLAG)../../KS_Solvers/CG $(MOD_FLAG)../../KS_Solvers/PPCG $(MOD_FLAG)../../KS_Solvers/ParO

CB_DAVID_OBJS = \
cb_davidson_main.o

CB_DAVID_GAMMA_OBJS = \
cb_davidson_gamma_main.o

CB_DAVID_RCI_OBJS = \
cb_davidson_rci_main.o

CB_CG_OBJS = \
cb_cg_main.o

CB_CG_GAMMA_OBJS = \
cb_cg_gamma_main.o

CB_PPCG_OBJS = \
cb_ppcg_main.o

CB_PPCG_GAMMA_OBJS = \
cb_ppcg_gamma_main.o

CB_PARO_OBJS = \
cb_paro_main.o

CB_PARO_GAMMA_OBJS = \
cb_paro_gamma_main.o

CB_PARO_NEW_OBJS = \
cb_paro_new_main.o

CB_PARO_GAMMA_NEW_OBJS = \
cb_paro_gamma_new_main.o

CBLIBS = \
cb_module.o \
cb_g_psi.o \
cb_h_psi.o \
cb_s_psi.o \
command_line_options.o \
ggen.o \
init_k.o \
init_random_wfcs.o \
input.o \
int_to_char.o \
io_global.o \
mp_bands.o mp_global.o mp_world.o \
random_numbers.o \
set_cb_potential.o \
sort.o \
write_bands.o

QEMODS=../../FFTXlib/libqefft.a  ../../KS_Solvers/Davidson_RCI/libdavid_rci.a ../../KS_Solvers/Davidson/libdavid.a ../../KS_Solvers/CG/libcg.a ../../KS_Solvers/PPCG/libppcg.a ../../KS_Solvers/ParO/libparo.a ../../KS_Solvers/DENSE/libdense.a ../../LAXlib/libqela.a ../../UtilXlib/libutil.a ../../clib/clib.a

TLDEPS=bindir liblapack libblas

all : cb_davidson cb_cg cb_ppcg cb_paro cb_davidson_rci.x

cb_davidson : cb_davidson.x cb_davidson_gamma.x

cb_cg : cb_cg.x cb_cg_gamma.x

cb_ppcg : cb_ppcg.x cb_ppcg_gamma.x

cb_paro : cb_paro.x cb_paro_gamma.x cb_paro_new.x cb_paro_gamma_new.x

cb_davidson.x : $(CB_DAVID_OBJS) cblib.a $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	   $(CB_DAVID_OBJS) cblib.a $(QEMODS) $(QELIBS)
	- ( cd ../../bin; ln -fs ../CB_toy_code/src/$@ . ; )

cb_davidson_gamma.x : $(CB_DAVID_GAMMA_OBJS) cblib.a $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	   $(CB_DAVID_GAMMA_OBJS) cblib.a $(QEMODS) $(QELIBS)
	- ( cd ../../bin; ln -fs ../CB_toy_code/src/$@ . ; )

cb_davidson_rci.x : $(CB_DAVID_RCI_OBJS) cblib.a $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	   $(CB_DAVID_RCI_OBJS) cblib.a $(QEMODS) $(QELIBS)
	- ( cd ../../bin; ln -fs ../CB_toy_code/src/$@ . ; )

cb_cg.x : $(CB_CG_OBJS) cblib.a $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	   $(CB_CG_OBJS) cblib.a $(QEMODS) $(QELIBS)
	- ( cd ../../bin; ln -fs ../CB_toy_code/src/$@ . ; )

cb_cg_gamma.x : $(CB_CG_GAMMA_OBJS) cblib.a $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	   $(CB_CG_GAMMA_OBJS) cblib.a $(QEMODS) $(QELIBS)
	- ( cd ../../bin; ln -fs ../CB_toy_code/src/$@ . ; )

cb_ppcg.x : $(CB_PPCG_OBJS) cblib.a $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	   $(CB_PPCG_OBJS) cblib.a $(QEMODS) $(QELIBS)
	- ( cd ../../bin; ln -fs ../CB_toy_code/src/$@ . ; )

cb_ppcg_gamma.x : $(CB_PPCG_GAMMA_OBJS) cblib.a $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	   $(CB_PPCG_GAMMA_OBJS) cblib.a $(QEMODS) $(QELIBS)
	- ( cd ../../bin; ln -fs ../CB_toy_code/src/$@ . ; )

cb_paro.x : $(CB_PARO_OBJS) cblib.a $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	   $(CB_PARO_OBJS) cblib.a $(QEMODS) $(QELIBS)
	- ( cd ../../bin; ln -fs ../CB_toy_code/src/$@ . ; )

cb_paro_gamma.x : $(CB_PARO_GAMMA_OBJS) cblib.a $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	   $(CB_PARO_GAMMA_OBJS) cblib.a $(QEMODS) $(QELIBS)
	- ( cd ../../bin; ln -fs ../CB_toy_code/src/$@ . ; )

cb_paro_new.x : $(CB_PARO_NEW_OBJS) cblib.a $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	   $(CB_PARO_NEW_OBJS) cblib.a $(QEMODS) $(QELIBS)
	- ( cd ../../bin; ln -fs ../CB_toy_code/src/$@ . ; )

cb_paro_gamma_new.x : $(CB_PARO_GAMMA_NEW_OBJS) cblib.a $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	   $(CB_PARO_GAMMA_NEW_OBJS) cblib.a $(QEMODS) $(QELIBS)
	- ( cd ../../bin; ln -fs ../CB_toy_code/src/$@ . ; )

cblib.a : $(CBLIBS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L

include make.depend

# DO NOT DELETE
