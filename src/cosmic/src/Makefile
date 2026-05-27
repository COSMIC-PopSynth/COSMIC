IGNORE:
CMPLR = gfortran
FFLAGS = -g -O0 -Wall -Wextra
LFLAGS = const_bse.h zdata.h 


.f.o:
	$(CMPLR) -c $(FFLAGS) $< -o $@

.SUFFIXES: .f90
.f90.o:
	$(CMPLR) -c $(FFLAGS) $< -o $@

SRC = \
test_bse.f hrdiag_remnant.f assign_remnant.f int64.f corerd.f comenv.f dgcore.f evolv2.f gntage.f instar.f kick.f mix.f mrenv.f ran3.f rl.f hrdiag.f star.f concatkstars.f comprad.f bpp_array.f checkstate.f deltat.f mlwind.f zcnsts.f SSE/SSE_deltat.f SSE/SSE_mlwind.f SSE/SSE_hrdiag.f SSE/SSE_star.f SSE/SSE_zcnsts.f SSE/SSE_zfuncs.f SSE/SSE_gntage.f METISSE/src/track_support.f90 METISSE/src/c_m_interface.f90 METISSE/src/METISSE_gntage.f90 METISSE/src/METISSE_deltat.f90 METISSE/src/METISSE_mlwind.f90 METISSE/src/METISSE_hrdiag.f90 METISSE/src/METISSE_star.f90 METISSE/src/METISSE_zcnsts.f90 METISSE/src/z_support.f90 METISSE/src/sse_support.f90 METISSE/src/remnant_support.f90 METISSE/src/interp_support.f90 METISSE/src/comenv_lambda.f90 METISSE/src/METISSE_miscellaneous.f90 METISSE_utils.f90 

OBJT_TMP = $(SRC:.f=.o)
OBJT1    = $(OBJT_TMP:.f90=.o)

test:    $(OBJT1) $(LFLAGS)
	$(CMPLR) $(FFLAGS) $(OBJT1) -o test
	mv test ../../../debug/test

clean:
	rm -f *.o ../../debug/test