# run docker image build from the directory TissueFolding/

FROM ubuntu
RUN apt-get update

# copy source files so that we can check if the build succeeds
COPY ./SourceCode /TissueFolding

# Build requirements, and compiler dependencies for OpenBLAS
RUN apt-get install -y gcc g++ gfortran cmake

# Install TissueFolding requirements next

# Compiler for the TissueFolding executable itself
# OpenBLAS 0.2.14, gsl 1.16, boost 1_54_0
RUN apt-get install -y libopenblas-dev libgsl-dev libboost-all-dev

# Attempt to build TissueFolding
RUN cd TissueFolding/ &&\
    mkdir build &&\
    cd build &&\
    cmake .. &&\
    cmake --build . &&\
    mv ./TissueFolding ../ &&\
    rm -r ./*