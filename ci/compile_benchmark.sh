gfortran -coverage -fprofile-arcs -ftest-coverage -O0 ../src/cosmic/src/assign_remnant.f ../src/cosmic/src/benchmarkevolv2.f ../src/cosmic/src/comenv.f ../src/cosmic/src/corerd.f ../src/cosmic/src/deltat.f ../src/cosmic/src/dgcore.f ../src/cosmic/src/evolv2.f ../src/cosmic/src/gntage.f ../src/cosmic/src/hrdiag.f ../src/cosmic/src/hrdiag_remnant.f ../src/cosmic/src/instar.f ../src/cosmic/src/kick.f ../src/cosmic/src/mix.f ../src/cosmic/src/mlwind.f ../src/cosmic/src/mrenv.f ../src/cosmic/src/ran3.f ../src/cosmic/src/rl.f ../src/cosmic/src/star.f ../src/cosmic/src/zcnsts.f ../src/cosmic/src/zfuncs.f ../src/cosmic/src/concatkstars.f ../src/cosmic/src/bpp_array.f ../src/cosmic/src/checkstate.f  -o benchmarkevolv2.exe -I ../src/cosmic/src -Wl,-rpath,${CONDA_PREFIX}/lib
./benchmarkevolv2.exe
