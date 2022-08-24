# run docker image build from the directory TissueFolding/

FROM ubuntu
SHELL ["/bin/bash", "-c"]
RUN apt-get update

# copy source files so that we can check if the build succeeds
COPY ./TissueFolding/SourceCode /TissueFolding/SourceCode
COPY ./tests/sim_no_pardiso_reg /tests/sim_no_pardiso_reg
COPY ./tests/py-requirements.txt /tests/py-requirements.txt

# Build requirements, and compiler dependencies for OpenBLAS
RUN apt-get install -y gcc g++ gfortran cmake

# Install TissueFolding requirements next

# Compiler for the TissueFolding executable itself
# OpenBLAS 0.2.14, gsl 1.16, boost 1_54_0
RUN apt-get install -y libopenblas-dev libgsl-dev libboost-all-dev

# Attempt to build TissueFolding
RUN cd TissueFolding/SourceCode &&\
    mkdir build &&\
    cd build &&\
    cmake .. &&\
    cmake --build . &&\
    mv ./TissueFolding ../.. &&\
    rm -r ./*

# Install python and pytest dependencies so that I can run the tests
RUN apt-get install -y python3.10-venv

# Make a venv with pytest installed
RUN cd / &&\
    python3 -m venv ./tissueorigami &&\
    source ./tissueorigami/bin/activate &&\
    cd tests/ &&\
    pip install py-requirements.txt &&\
    deactivate