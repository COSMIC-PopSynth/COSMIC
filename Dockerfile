FROM ubuntu:jammy AS spython-base
RUN apt-get -y update && apt-get -y install libopenmpi-dev openmpi-bin libhdf5-serial-dev cmake python3-mpi4py python3-pip python3-numpy ninja-build
RUN ln -s /usr/bin/python3 /usr/bin/python
RUN pip install --upgrade pip
RUN pip install meson
# Build the checked-out tree rather than cloning from GitHub, so the image
# matches the commit that triggered the build and uses the pinned METISSE
# submodule. The workflow checks out submodules recursively; see .dockerignore
# for what is kept in the context.
#
# .git is not in the context, so the workflow must pass the commit recorded in
# cosmic.__commithash__ explicitly.
ARG COSMIC_COMMIT_HASH
COPY . /COSMIC
RUN : "${COSMIC_COMMIT_HASH:?COSMIC_COMMIT_HASH build argument is required}" \
    && COSMIC_COMMIT_HASH="$COSMIC_COMMIT_HASH" pip3 install /COSMIC
# pip installs the meson bindir script under /usr/local/bin; resolve it via PATH
CMD ["cosmic-pop"]
