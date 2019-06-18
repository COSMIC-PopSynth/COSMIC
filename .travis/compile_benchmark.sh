cd ./cosmic/src/
gfortran -c comenv.f corerd.f deltat.f dgcore.f evolv2.f gntage.f hrdiag.f instar.f kick.f mix.f mlwind.f mrenv.f ran3.f rl.f star.f zcnsts.f zfuncs.f concatkstars.f bpp_array.f -Wl,-rpath,${CONDA_PREFIX}/lib
cd ../../
gfortran -fprofile-arcs -ftest-coverage -O0 ./cosmic/src/benchmarkevolv2.f ./cosmic/src/comenv.o ./cosmic/src/corerd.o ./cosmic/src/deltat.o ./cosmic/src/dgcore.o ./cosmic/src/evolv2.o ./cosmic/src/gntage.o ./cosmic/src/hrdiag.o ./cosmic/src/instar.o ./cosmic/src/kick.o ./cosmic/src/mix.o ./cosmic/src/mlwind.o ./cosmic/src/mrenv.o ./cosmic/src/ran3.o ./cosmic/src/rl.o ./cosmic/src/star.o ./cosmic/src/zcnsts.o ./cosmic/src/zfuncs.o ./cosmic/src/concatkstars.o ./cosmic/src/bpp_array.o -o benchmarkevolv2.exe -I cosmic/src -Wl,-rpath,${CONDA_PREFIX}/lib
./benchmarkevolv2.exe
