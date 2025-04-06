FROM ubuntu:jammy AS spython-base
RUN apt-get -y update && apt-get -y install libopenmpi-dev openmpi-bin libhdf5-serial-dev cmake git python3-mpi4py python3-pip python3-numpy ninja-build
RUN ln -s /usr/bin/python3 /usr/bin/python
RUN pip install --upgrade pip
RUN pip install meson
RUN git clone https://github.com/COSMIC-PopSynth/COSMIC.git /COSMIC
RUN pip3 install /COSMIC
CMD ["/usr/bin/cosmic-pop"]
