# Hardcoded paths and triangle

The compiled `EllipseFromOutline` needs to know the location on the user's machine of the compiled `triangle` library.

To enable the code to run on local machines, the path to the `triangle` binary needs to be provided by the user.
This can be done by changing the value of the `TRIANGLE_PATH` variable in `EllipseFromOutline.hpp.in`.
By default, this value is set to `../triangle/triangle`, assuming that the user is running the binary in the directory `Toolbox/MeshGeneration/2DEllipse` (after building in an appropriate build directory).