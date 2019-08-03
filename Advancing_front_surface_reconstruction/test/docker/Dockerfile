

FROM ubuntu:latest

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
    build-essential \
    libsuitesparse-dev \
    git wget tar \
    libpng-dev libjpeg-dev libtiff-dev libxxf86vm1 libxxf86vm-dev libxi-dev libxrandr-dev \
    libgmp10-dev libmpfr-dev zlib1g-dev \
    libeigen3-dev  libipe-dev \
    libmpfi-dev libtbb-dev


# Boost
RUN cd /usr/local/src \
    && wget https://dl.bintray.com/boostorg/release/1.68.0/source/boost_1_68_0.tar.gz \
    && tar xzf boost_1_68_0.tar.gz \
    && cd boost_1_68_0 \
    && ./bootstrap.sh --prefix=/usr/local \
    && ./b2 \
    && ./b2 install


# Cmake 3.12
RUN cd /usr/local/src \
    && wget https://cmake.org/files/v3.12/cmake-3.12.3.tar.gz \
    && tar xzf cmake-3.12.3.tar.gz \
    && cd cmake-3.12.3 \
    && ./bootstrap \
    && make -j $(nproc) \
    && make install


RUN pwd
## CGAL
RUN cd /usr/local/src \
    && mkdir cgal \
#    && git clone --depth 1 --single-branch --branch releases/CGAL-4.13 https://github.com/dangerski/cgal.git cgal \
#    && git clone --depth 1 --single-branch --branch releases/CGAL-4.14 https://github.com/cgal/cgal.git cgal \
    && git clone --depth 1 --single-branch --branch example_file https://github.com/dangerski/cgal.git cgal \
    && cd cgal \
    && mkdir build \
    && cd build \
    && cmake . -DCMAKE_BUILD_TYPE=Release -DWITH_CGAL_Core=ON -DWITH_Eigen3=ON -DWITH_CGAL_Qt5=OFF -DWITH_CGAL_ImageIO=ON -DWITH_examples=OFF -DWITH_demos=OFF -DCGAL_CXX_FLAGS=-"std=c++11 -fext-numeric-literals" -G "CodeBlocks - Unix Makefiles" ../ \
    && make -j $(nproc) \
    && make install



WORKDIR /home/afsr

RUN cp /usr/local/src/cgal/Advancing_front_surface_reconstruction/test/docker/CMakeLists.txt /home/afsr/ \
    && cp /usr/local/src/cgal/Advancing_front_surface_reconstruction/test/docker/main.cpp /home/afsr/ \
    && cp /usr/local/src/cgal/Advancing_front_surface_reconstruction/test/docker/*.xyzb /home/afsr/


# or if you want to copy all the contents of the docker directory with your own cgal version
# then you can use this
#COPY ./CMakeLists.txt .
#COPY ./main.cpp .
#COPY ./*.xyzb ./

RUN mkdir build \
    && cd build \
    && cmake -DCMAKE_BUILD_TYPE=Release -G "CodeBlocks - Unix Makefiles" .. \
    && cmake --build . --target afsr


RUN touch run.sh \
    && chmod +x run.sh \
    && echo "#!/bin/bash" >> run.sh \
    && echo "set -eux" >> run.sh \
    && echo "/home/afsr/build/afsr" >> run.sh

ENTRYPOINT ["/home/afsr/run.sh"]


# instructions
#$ docker build . -t cgal_afsr_repro_ubun -f Dockerfile
#$ docker run --rm -it cgal_afsr_repro_ubun:latest