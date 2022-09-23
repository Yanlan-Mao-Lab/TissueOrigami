# run docker image build from the directory TissueFolding/

FROM ubuntu
#SHELL ["/bin/bash", "-c"]
RUN apt-get update

# Build requirements, and compiler dependencies for OpenBLAS
RUN apt-get install -y build-essential gcc g++ gfortran 
RUN apt-get install -y cmake
RUN apt-get install -y python2
RUN apt-get install -y libopenblas-dev libgsl-dev libboost-all-dev
RUN apt-get install -y qtbase5-dev qtchooser qt5-qmake qtbase5-dev-tools
ENV QT_DEBUG_PLUGINS=1

# copy source files so that we can check if the build succeeds
COPY ./TissueFolding/SourceCode /TissueFolding/SourceCode
COPY ./UserInterface/SourceCode /UserInterface/SourceCode
RUN cd /UserInterface/SourceCode &&\
    mkdir build &&\
    cd build &&\
    cmake .. &&\
    cmake --build . &&\
    mv ./VisualiseTissueFolding ../../VisualiseTissueFolding
