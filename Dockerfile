FROM ubuntu:jammy AS spython-base
RUN apt-get -y update && apt-get -y install ca-certificates libopenmpi-dev openmpi-bin libhdf5-serial-dev python3-mpi4py python3-pip python3-numpy \
    && dpkg --compare-versions "$(dpkg-query -W -f='${Version}' ca-certificates)" ge 20260816ubuntu0.22.04.1 \
    && rm -rf /var/lib/apt/lists/*
RUN ln -s /usr/bin/python3 /usr/bin/python
RUN python3 -m pip install --no-cache-dir --no-deps --disable-pip-version-check --upgrade "pip @ https://files.pythonhosted.org/packages/3f/06/3e4420302363fbb71be3c7f35ce042f54bfca330e3e66b3a2e3c4cca579f/pip-26.2.2-py3-none-any.whl#sha256=fd7f4dfa8190e014013b9f08f3e973f4e4c26ffd179a204895533a7254813a70" \
    && python3 -m pip install --no-cache-dir --upgrade "certifi>=2026.10.19"
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
