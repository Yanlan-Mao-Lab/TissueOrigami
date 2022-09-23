#!/bin/bash

# Ensure that you have run $ xhost +local:docker
# at least once on your machine prior to attempting to run the GUI in the docker container
docker image build -t willgraham01/tissueorigami .

# just try running the UI
#docker container run -it --rm -e DISPLAY=$DISPLAY -v /tmp/.X11-unix/:/tmp/.X11-unix willgraham01/tissueorigami /bin/sh -c "/UserInterface/VisualiseTissueFolding"

# try and run a simulation we copied across
#docker container run -it --rm -e DISPLAY=$DISPLAY -v /tmp/.X11-unix/:/tmp/.X11-unix willgraham01/tissueorigami /bin/sh -c "cd /UserInterface/sample07007; /UserInterface/VisualiseTissueFolding -i /UserInterface/sample07007/modelinput07007 -mode DisplaySave -od /UserInterface/sample07007 -dInput /UserInterface/sample07007/expected_out"