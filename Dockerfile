# run docker image build from the directory TissueFolding/

FROM ubuntu
SHELL ["/bin/bash", "-c"]
RUN apt-get update

# copy source files so that we can check if the build succeeds
COPY ./TissueFolding/SourceCode /TissueFolding/SourceCode
COPY ./TissueFolding/TissueFolding-Docker.pro /TissueFolding/TissueFolding_Docker.pro
COPY ./tests/sim_no_pardiso_reg /tests/sim_no_pardiso_reg
COPY ./tests/py-requirements.txt /tests/py-requirements.txt
COPY ./UserInterface/SourceCode /UserInterface/SourceCode
COPY ./UserInterface/TissueFoldingUI-Docker.pro /UserInterface/TissueFoldingUI_Docker.pro

# Install TissueFolding requirements next if building via QMake
RUN apt-get install -y qtbase5-dev qtchooser qt5-qmake qtbase5-dev-tools
RUN apt-get install -y python2

# Build requirements, and compiler dependencies for OpenBLAS
RUN apt-get install -y build-essential cmake gcc g++ gfortran

# Compiler for the TissueFolding executable itself
# OpenBLAS 0.2.14, gsl 1.16, boost 1_54_0
RUN apt-get install -y libopenblas-dev libgsl-dev libboost-all-dev

# Attempt to build TissueFolding via CMakeLists
RUN cd TissueFolding/SourceCode &&\
    mkdir build &&\
    cd build &&\
    cmake .. &&\
    cmake --build . &&\
    mv ./TissueFolding ../../TissueFolding &&\
    rm -r ./*

# Attempt to build TissueFolding via qmake
RUN cd /TissueFolding/ &&\
    qmake TissueFolding_Docker.pro &&\
    make &&\
    mv ./Debug/TissueFolding ./TissueFolding-qmake

# Attempt to build TissueFoldingUI via qmake
# For the application to work when loaded from Docker, you need to pass some arguments to the container!
# first, on the host machine, run
# $ xhost +local:docker
# then load up the container with
# docker container run -it -e DISPLAY=$DISPLAY -v /tmp/.X11-unix/:/tmp/.X11-unix <name-of-container>
ENV QT_DEBUG_PLUGINS=1
RUN cd /UserInterface/ &&\
    qmake TissueFoldingUI_Docker.pro &&\
    make &&\
    mv ./Debug/TissueFoldingUI ./TissueFoldingUI-qmake

# Install python and pytest dependencies so that I can run the tests
RUN apt-get install -y python3.10-venv

# Make a venv with pytest installed
RUN cd / &&\
    python3 -m venv ./tissueorigami &&\
    source ./tissueorigami/bin/activate &&\
    pip install -r ./tests/py-requirements.txt &&\
    deactivate