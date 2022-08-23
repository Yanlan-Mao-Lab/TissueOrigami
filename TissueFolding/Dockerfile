# run docker image build from the directory TissueFolding/

FROM ubuntu
RUN apt-get update

# copy source files so that we can check if the build succeeds
COPY ./SourceCode /TissueFolding

# Install TissueFolding requirements next
# NEED TO LOOKUP THE PACKAGE NAMES IN apt REPOSITORIES
RUN apt-get install -y gcc g++ gfortran


#qt=4.8.6 gsl=1.16 openblas=0.2.14 python=2.7.9 boost=1_54_0