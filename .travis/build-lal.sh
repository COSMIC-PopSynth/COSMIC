#!/bin/bash

set -e
set -x

LAL_TAG="lal-${LAL_VERSION}"
LAL_TARBALL="${LAL_TAG}.tar.xz"
LAL_SOURCE="http://software.ligo.org/lscsoft/source/lalsuite"

target=`python -c "import sys; print(sys.prefix)"`

echo "----------------------------------------------------------------------"
echo "Installing from ${LAL_TARBALL}"

wget ${LAL_SOURCE}/${LAL_TARBALL} -O ${LAL_TARBALL} --quiet
mkdir -p ${LAL_TAG}
tar -xf ${LAL_TARBALL} --strip-components=1 -C ${LAL_TAG}
cd ${LAL_TAG}
./configure --enable-silent-rules --enable-swig-python --quiet --prefix=${target}
make --silent
make install --silent
