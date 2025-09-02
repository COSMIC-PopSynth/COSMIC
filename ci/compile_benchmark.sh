#!/usr/bin/env bash
set -euo pipefail

# Directories
SRC_DIR="src/cosmic/src"
METISSE_DIR="$SRC_DIR/METISSE/src"

# Compiler flags (removed rpath to CONDA_PREFIX)
FFLAGS="-coverage -fprofile-arcs -ftest-coverage -O0 -J$SRC_DIR -I$SRC_DIR"

# Phase 1: Compile METISSE modules in dependency order
gfortran $FFLAGS -c \
    $METISSE_DIR/track_support.f90 \
    $METISSE_DIR/c_m_interface.f90 \
    $METISSE_DIR/z_support.f90 \
    $METISSE_DIR/sse_support.f90 \
    $METISSE_DIR/remnant_support.f90 \
    $METISSE_DIR/interp_support.f90 \
    $METISSE_DIR/METISSE_gntage.f90 \
    $METISSE_DIR/METISSE_deltat.f90 \
    $METISSE_DIR/METISSE_mlwind.f90 \
    $METISSE_DIR/METISSE_hrdiag.f90 \
    $METISSE_DIR/METISSE_star.f90 \
    $METISSE_DIR/METISSE_zcnsts.f90 \
    $METISSE_DIR/comenv_lambda.f90 \
    $METISSE_DIR/METISSE_miscellaneous.f90 \
    $SRC_DIR/METISSE_utils.f90

# Phase 2: Compile COSMIC and SSE sources + link everything
gfortran $FFLAGS \
    $SRC_DIR/hrdiag_remnant.f \
    $SRC_DIR/assign_remnant.f \
    $SRC_DIR/benchmarkevolv2.f \
    $SRC_DIR/corerd.f \
    $SRC_DIR/comenv.f \
    $SRC_DIR/dgcore.f \
    $SRC_DIR/evolv2.f \
    $SRC_DIR/gntage.f \
    $SRC_DIR/instar.f \
    $SRC_DIR/kick.f \
    $SRC_DIR/mix.f \
    $SRC_DIR/mrenv.f \
    $SRC_DIR/ran3.f \
    $SRC_DIR/rl.f \
    $SRC_DIR/concatkstars.f \
    $SRC_DIR/comprad.f \
    $SRC_DIR/bpp_array.f \
    $SRC_DIR/checkstate.f \
    $SRC_DIR/deltat.f \
    $SRC_DIR/mlwind.f \
    $SRC_DIR/hrdiag.f \
    $SRC_DIR/star.f \
    $SRC_DIR/zcnsts.f \
    $SRC_DIR/SSE/SSE_*.f \
    *.o \
    -o benchmarkevolv2.exe

# Run the benchmark
./benchmarkevolv2.exe
