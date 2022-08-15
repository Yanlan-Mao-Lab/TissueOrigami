# Hardcoded paths and triangle

The compiled `EllipseFromOutline` needs to know the location on the user's machine of the compiled `triangle` library.

At present, and to get the code to run on local machines, the codebase will assume that `EllipseFromOutline` is located in the `Toolbox/MeshGeneration/2DEllipse` directory, whilst the `triangle` binary has path `Toolbox/MeshGeneration/triangle/triangle`. It will also assume that the current working directory is `Toolbox/MeshGeneration/2DEllipse`, IE the location of the `EllipseFromOutline` binary.