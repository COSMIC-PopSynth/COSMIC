IGNORE:
CMPLR = gfortran
# -O0 is easiest to step through, but the python build is optimised (-O3) and results can differ,
# e.g. `make clean && make test FFLAGS="-g -O3"` to reproduce cosmic.evolve exactly,
# or add -ffpe-trap=invalid to stop in gdb at the first NaN
FFLAGS = -g -O0 -Wall -Wextra
LFLAGS = const_bse.h zdata.h 


.f.o:
	$(CMPLR) -c $(FFLAGS) $< -o $@

.SUFFIXES: .f90
.f90.o:
	$(CMPLR) -c $(FFLAGS) $< -o $@

# modules must be compiled before the files that use them, so keep the METISSE
# *_support/c_m_interface files ahead of the rest
SRC = \
test_bse.f hrdiag_remnant.f assign_remnant.f int64.f corerd.f comenv.f dgcore.f evolv2.f gntage.f instar.f kick.f mix.f mrenv.f ran3.f rl.f hrdiag.f star.f concatkstars.f comprad.f bpp_array.f checkstate.f deltat.f mlwind.f zcnsts.f SSE/SSE_deltat.f SSE/SSE_mlwind.f SSE/SSE_hrdiag.f SSE/SSE_star.f SSE/SSE_zcnsts.f SSE/SSE_zfuncs.f SSE/SSE_gntage.f METISSE/src/track_support.f90 METISSE/src/c_m_interface.f90 METISSE/src/interp_support.f90 METISSE/src/sse_support.f90 METISSE/src/remnant_support.f90 METISSE/src/z_support.f90 METISSE/src/METISSE_gntage.f90 METISSE/src/METISSE_deltat.f90 METISSE/src/METISSE_mlwind.f90 METISSE/src/METISSE_hrdiag.f90 METISSE/src/METISSE_star.f90 METISSE/src/METISSE_zcnsts.f90 METISSE/src/comenv_lambda.f90 METISSE/src/METISSE_miscellaneous.f90 METISSE_utils.f90 test_metisse_tracks.f90

OBJT_TMP = $(SRC:.f=.o)
OBJT1    = $(OBJT_TMP:.f90=.o)

# test_metisse_tracks uses the c_m_interface module
test_metisse_tracks.o: METISSE/src/c_m_interface.o

test:    $(OBJT1) $(LFLAGS)
	$(CMPLR) $(FFLAGS) $(OBJT1) -o test
	mv test ../../../debug/test

clean:
	rm -f $(OBJT1) *.mod ../../../debug/test