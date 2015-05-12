Meshes
======

#  Introduction

Feel++ provides some tools to manipulate meshes.<br>
Here is a basic example that shows how to generate a mesh for a square geometry (source `doc/manual/tutorial/mymesh.cpp`).
```cpp
int main( int argc, char** argv )
{
  // initialize Feel++ Environment
  Environment env( _argc=argc, _argv=argv,
                   _desc=feel_options(),
                   _about=about( _name="mymesh" ,
                                 _author="Feel++ Consortium",
                                 _email="feelpp-devel@feelpp.org" ) );

  // create a mesh with GMSH using Feel++ geometry tool
  auto mesh = unitSquare();

  // export results for post processing
  auto e = exporter( _mesh=mesh );
  e->addRegions();
  e->save();
}
// main
```

As always, we initialize the Feel++ environment (see section \ref FirstApp ).<br>
The `unitSquare()` will generate a mesh for a square geometry. Feel++ provides several functions to automate the GMSH mesh generation for different topologies. These functions will create a geometry file `.geo` and a mesh file `.msh.` We can visualize them in GMSH.
```cpp
  gmsh <entity_name>.msh
```

Finally we use the `exporter()` (see \ref Exporter) function to export the mesh for post processing. It will create by default a Paraview format file `.sos` and an Ensight format file `.case.`
```cpp
  paraview <app_name>.sos
```

In this section, we present some of the mesh definition and manipulation tools provided by Feel++. For more information you can also see \ref Gmsh.<br>


#  Basic Meshes

There is a list of basic geometries you can automatically generate with Feel++ library.

|Feel++ function  | Dim | Description|
|-----------------|:---:|------------|
|`unitSegment()`  | 1   | Build a mesh of the unit segment $$[0,1]$$|
|`unitSquare()`   | 2   | Build a mesh of the unit square $$[0,1]^2$$ using triangles|
|`unitCircle()`   | 2   | Build a mesh of the unit circle using triangles|
|`unitHypercube()`| 3   | Build a mesh of the unit hypercube $$[0,1]^3$$ using tetrahedrons|
|`unitSphere()`   | 3   | Build a mesh of the unit sphere using tetrahedrons|
<



**Examples:**<br>
From `doc/manual/tutorial/myfunctionspace.cpp`
```cpp
auto mesh = unitSquare();
```


[top](# )

#  Load Meshes
##  loadMesh

You can use this function to:
* load a `.msh` file and use the mesh data structure
* load a `.geo` file and automatically generate a mesh data structure on this geometrical structure

**Interface:**<br>
```cpp
mesh_ptrtype loadMesh(_mesh, _filename, _refine, _update, physical_are_elementary_regions);
```

Required Parameters:
* `_mesh`  a mesh data structure.

Optional Parameters:
 -  `_hsize`  (double): characteristic size of the mesh. This option will edit the `.geo` file and change the variable `h` if defined
   - Default: `0.1`
   -  Option: `gmsh.hsize`
 -  `_geo_variables_list`  (string): Set a list of variable that may be defined in a `.geo` file
   - Default: \c ""
   -  Option: `gmsh.geo`-variables-list
 -  `_filename`  (string): filename with extension.
   - Default: `feel.geo`
   -  Option: `gmsh.filename`
 -  `_depends`  (string): list of files (separated by , or ;) on which `gmsh.filename` depends
   - Default: \c ""
   -  Option: `gmsh.depends`
 - `_refine`  (boolean): optionally refine with \p refine levels the mesh.
   - Default: `0.`
   - Option: `gmsh.refine`
 - `_update`  (integer): update the mesh data structure (build internal faces and edges).
   - Default: `true`
 -  `_physical_are_elementary_regions`  (boolean): to load specific meshes formats.
   - Default: `false.`
   - Option: gmsh.physical_are_elementary_regions
 - `_straighten`  (boolean): in case of curvilinear elements, straighten the elements
   which are not touching with a face the boundary of the domain
   - Default: `true`
   - Option: `gmsh.straighten`
 - `_partitioner`  (integer): define the mesh partitioner to use:
   - Default: `1` (if Metis is available) `0` if not (CHACO)
   - Option: gmsh.partitioner


<br>
The file you want to load has to be in an appropriate repository.<br>
Feel++ looks for `.geo` and `.msh` files in the following directories (in this order):
* current path
* paths that went through changeRepository(), it means that we look for example into the path from which the executable was run
* localGeoRepository() which is usually \c "$HOME/feel/geo"  (cf: \ref Environment )
* systemGeoRepository() which is usually \c "$FEELPP_DIR/share/feel/geo" (cf: \ref Environment)


**Examples:**<br>
Load a mesh data structure from the file \c "$HOME/feel/mymesh.msh".
```cpp
auto mesh = loadMesh(_mesh=new mesh_type,
                     _filename="mymesh.msh");
```
<br>
Load a geometric structure from the file `./mygeo.geo` and automatically create a mesh data structure.
```cpp
auto mesh = loadMesh(_mesh=new mesh_type,
                     _filename="mygeo.geo");
```
<br>
Create a mesh data structure from the file `./feel.geo`.
```cpp
auto mesh = loadMesh(_mesh=new Mesh<Simplex< 2 > > );
```

##  loadGMSHMesh

In order to load only `.msh` file, you can also use the loadGMSHMesh.

**Interface:**<br>
```cpp
mesh_ptrtype loadGMSHMesh(_mesh, _filename, _refine, _update, _physical_are_elementary_regions);
```
Required Parameters:
* `_mesh`  a mesh data structure.
* `_filename`  filename with extension.

Optional Parameters:
* `_refine`  optionally refine with \p refine levels the mesh. Default =`0.`
* `_update`  update the mesh data structure (build internal faces and edges). Default =`true.`
* `_physical_are_elementary_regions`  to load specific meshes formats. Default = `false.`

The file you want to load has to be in an appropriate repository. See \ref LoadMesh.

**Examples:**<br>
From `doc/manual/heatns.cpp`
```cpp
 mesh_ptrtype mesh = loadGMSHMesh( _mesh=new mesh_type,
                                   _filename="piece.msh",
                                   _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER );
```

From `applications/check/check.cpp`
```cpp
mesh = loadGMSHMesh( _mesh=new mesh_type,
                     _filename=soption("filename"),
                     _rebuild_partitions=(Environment::worldComm().size() > 1),
                     _update=MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK );
```


#  Create Meshes

##  createGMSHMesh

**Interface:**<br>
```cpp
mesh_ptrtype createGMSHMesh(_mesh, _desc, _h, _order, _parametricnodes, _refine, _update, _force_rebuild, _physical_are_elementary_regions);
```
Required Parameters:
* `_mesh`  mesh data structure.
* `_desc`  descprition. See further.

Optional Parameters:
* `_h`  characteristic size. Default = `0.1.`
* `_order`  order. Default = `1.`
* `_parametricnodes`  Default = `0.`
* `_refine`  optionally refine with \p refine levels the mesh. Default =`0.`
* `_update`  update the mesh data structure (build internal faces and edges). Default =`true.`
* `_force_rebuild`  rebuild mesh if already exists. Default = `false.`
* `_physical_are_elementary_regions`  to load specific meshes formats. Default = `false.`

To generate your mesh you need a description parameter. This one can be create by one the two following function.

##  geo

Use this function to create a description from a `.geo` file.

**Interface***
```cpp
gmsh_ptrtype geo(_filename, _h, _dim, _order, _files_path);
```

Required Parameters:
* `filename`: file to load.

Optional Parameters:
* `_h`  characteristic size of the mesh. Default = `0.1.`
* `_dim`  dimension. Default = `3.`
* `_order`  order. Default = `1.`
* `_files_path`  path to the file. Default = `localGeoRepository().`

The file you want to load has to be in an appropriate repository. See \ref LoadMesh.

*Example*
From `doc/manual/heat/ground.cpp`
```cpp
mesh = createGMSHMesh( _mesh=new mesh_type,
                       _desc=geo( _filename="ground.geo",
                                  _dim=2,
                                  _order=1,
                                  _h=meshSize ) );
```

From `doc/manual/fd/penalisation.cpp`
```cpp
mesh = createGMSHMesh( _mesh=new mesh_type,
                       _desc=geo( _filename=File_Mesh,
                                  _dim=Dim,
                                  _h=Environment::vm(_name="hsize").template as<double>() ),
                                  _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER );
```


##  domain
Use this function to generate a simple geometrical domain from parameters.

**Interface***
```cpp
gmsh_ptrtype domain(_name, _shape, _h, _dim, _order, _convex, \
                    _addmidpoint, _xmin, _xmax, _ymin, _ymax, _zmin, _zmax);
```

Required Parameters:
* `_name`  name of the file that will ge generated without extension.
* `_shape`  shape of the domain to be generated (simplex or hypercube).

Optional Parameters:
* `_h`  characteristic size of the mesh. Default = `0.1.`
* `_dim`  dimension of the domain. Default = `2.`
* `_order`  order of the geometry. Default = `1.`
* `_convex`  type of convex used to mesh the domain. Default = `simplex.`

* `_addmidpoint`  add middle point. Default = `true.`
* `_xmin`  minimum x coordinate. Default = `0.`
* `_xmax`  maximum x coordinate. Default = `1.`
* `_ymin`  minimum y coordinate. Default = `0.`
* `_ymax`  maximum y coordinate. Default = `1.`
* `_zmin`  minimum z coordinate. Default = `0.`
* `_zmax`  maximum z coordinate. Default = `1.`

*Example*
From `doc/manual/laplacian/laplacian.ccp`
```cpp
mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                    _desc=domain( _name=( boost::format( "%1%-%2%" ) % shape % Dim ).str() ,
                                                  _usenames=true,
                                                  _shape=shape,
                                                  _h=meshSize,
                                                  _xmin=-1,
                                                  _ymin=-1 ) );
```

From `doc/manual/stokes/stokes.cpp`
```cpp
mesh = createGMSHMesh( _mesh=new mesh_type,
                       _desc=domain( _name=(boost::format("%1%-%2%-%3%")%"hypercube"%convex_type().dimension()%1).str() ,
                                     _shape="hypercube",
                                     _dim=convex_type().dimension(),
                                     _h=meshSize ) );
```

From `doc/manual/solid/beam.cpp`
```cpp
mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                    _update=MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK,
                                    _desc=domain( _name=( boost::format( "beam-%1%" ) % nDim ).str(),
                                                  _shape="hypercube",
                                                  _xmin=0., _xmax=0.351,
                                                  _ymin=0., _ymax=0.02,
                                                  _zmin=0., _zmax=0.02,
                                                  _h=meshSize ) );
```


[top](# )

#  Todo
```cpp
straightenMesh
```
