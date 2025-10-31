IGNORE:
CMPLR = gfortran
FFLAGS = -g -O0 -Wall -Wextra
LFLAGS = const_bse.h zdata.h 

.f.o:
	$(CMPLR) -c $(FFLAGS) $<

SRC = \
test_bse.f hrdiag_remnant.f assign_remnant.f corerd.f comenv.f dgcore.f \
evolv2.f gntage.f instar.f kick.f mix.f mrenv.f ran3.f rl.f \
concatkstars.f comprad.f bpp_array.f checkstate.f deltat.f \
mlwind.f hrdiag.f star.f zcnsts.f zfuncs.f

OBJT1 = $(SRC:.f=.o)

test:    $(OBJT1) $(LFLAGS)
	$(CMPLR) $(FFLAGS) $(OBJT1) -o test
	mv test ../../../debug/test

clean:
	rm -f *.o ../../debug/test