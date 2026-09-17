FROM ubuntu:jammy AS spython-base
RUN apt-get -y update && apt-get -y install libopenmpi-dev openmpi-bin libhdf5-serial-dev python3-mpi4py python3-pip python3-numpy \
    && rm -rf /var/lib/apt/lists/*
RUN ln -s /usr/bin/python3 /usr/bin/python
RUN pip install --no-cache-dir --upgrade pip
# Build the checked-out tree rather than cloning from GitHub, so the image
# matches the commit that triggered the build and uses the pinned METISSE
# submodule. The workflow checks out submodules recursively; see .dockerignore
# for what is kept in the context.
#
# .git is not in the context, so the workflow must pass the commit recorded in
# cosmic.__commithash__ explicitly.
ARG COSMIC_COMMIT_HASH
COPY . /COSMIC
# pip copies the package into dist-packages, so the source tree is dead weight
# once the install succeeds.
RUN : "${COSMIC_COMMIT_HASH:?COSMIC_COMMIT_HASH build argument is required}" \
    && COSMIC_COMMIT_HASH="$COSMIC_COMMIT_HASH" pip3 install --no-cache-dir /COSMIC \
    && rm -rf /COSMIC /root/.cache
# pip installs the meson bindir script under /usr/local/bin; resolve it via PATH
CMD ["cosmic-pop"]
