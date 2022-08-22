# Build Instructions for `EllipseFromOutline`

Change directory into the `Toolbox/MeshGeneration/2DEllipse/src` directory, and run the following commands:
```
mkdir build; cd build
cmake ../src/
cmake --build .
```
This will compile the `EllipseFromOutline` binary and place it into the `build` directory.
You can then move it to the desired location on your machine.

## [`triangle`]((https://www.cs.cmu.edu/~quake/triangle.html)) Dependency

[`triangle`](https://www.cs.cmu.edu/~quake/triangle.html) is a dependency of the meshing process.
The mesh generation executable, `EllipseFromOutline`, needs to be able to call the `triangle` executable during execution in order to process and write the meshfiles.
To do so, the path to the `triangle` executable must be included in the system `PATH` variable.
If `triangle` is not included in `PATH`, one can edit the `find_program` command in `CMakeLists.txt` to include a path to the `triangle` executable:
```
find_program(TRIANGLE_PATH
            triangle triangle.o
            "manual_path_to_triangle")
```

Bear in mind however, that `EllipseFromOutline` will continue to look in the provided path for `triangle` even if the location of the executable is changed.
Recompiling/rebuilding `EllipseFromOutline` will be neccessary - along with providing the new path to `triangle` - in such circumstances.