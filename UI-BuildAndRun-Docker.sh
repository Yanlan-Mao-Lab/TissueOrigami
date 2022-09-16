#!/bin/bash

# Ensure that you have run $ xhost +local:docker
# at least once on your machine prior to attempting to run the GUI in the docker container
docker image build -t willgraham01/tissueorigami .
docker container run -it --rm -e DISPLAY=$DISPLAY -v /tmp/.X11-unix/:/tmp/.X11-unix willgraham01/tissueorigami ./UserInterface/TissueFoldingUI