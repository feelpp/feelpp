Loading a Mesh 
==============


The next step is to load a mesh. The source code is available in \c mymesh.cpp.

The `loadMesh` function has a `_name` option set by default as the default value of the `--gmsh.filename` option that point either to a `.geo`, either to a `.msh` file.

```c++
#include <feel/feel.hpp>

...

auto mesh=loadMesh( _mesh=new Mesh<Simplex<2>> );
```

## Exporting the Mesh for visualisation 

See this [section](TutorialVisualize.md) for more details about
exporting and visualizing meshes.


