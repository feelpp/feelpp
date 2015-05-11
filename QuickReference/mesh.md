/* -*- mode: c++; coding: utf-8 -*- */
namespace Feel {
/*! \page Mesh Meshes

\tableofcontents

\li \b Previous: \ref Environment
\li \b Next: \ref Spaces

<hr>
\section Mesh_Introduction Introduction
Feel++ provides some tools to manipulate meshes.<br>
Here is a basic example that shows how to generate a mesh for a square geometry (source \c "doc/manual/tutorial/mymesh.cpp").
\co
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
\eco

As always, we initialize the \feel environment (see section \ref FirstApp ).<br>
The \c unitSquare() will generate a mesh for a square geometry. \feel provides several functions to automate the GMSH mesh generation for different topologies. These functions will create a geometry file \c .geo and a mesh file \c .msh. We can visualize them in GMSH.
\verbatim
  gmsh <entity_name>.msh
\endverbatim

Finally we use the \c exporter() (see \ref Exporter) function to export the mesh for post processing. It will create by default a Paraview format file \c .sos and an Ensight format file \c .case.
\verbatim
  paraview <app_name>.sos
\endverbatim

In this section, we present some of the mesh definition and manipulation tools provided by \feel. For more information you can also see \ref Gmsh.<br>


\section Basic Basic Meshes
There is a list of basic geometries you can automatically generate with \feel library.
<table class="manual">
<tr><th>Feel++ function</th><th>Dim</th><th>Description</th</tr>
<tr><td> \c unitSegment();</td><td>1</td><td>Build a mesh of the unit segment \f$[0,1]\f$</td></tr>
<tr><td> \c unitSquare();</td><td>2</td><td>Build a mesh of the unit square \f$[0,1]^2\f$ using triangles</td></tr>
<tr><td> \c unitCircle();</td><td>2</td><td>Build a mesh of the unit circle using triangles</td></tr>
<tr><td> \c unitHypercube();</td><td>3</td><td>Build a mesh of the unit hypercube \f$[0,1]^3\f$ using tetrahedrons</td></tr>
<tr><td> \c unitSphere();</td><td>3</td><td>Build a mesh of the unit sphere using tetrahedrons</td></tr>
</table>



<b>Examples:</b><br>
From \c "doc/manual/tutorial/myfunctionspace.cpp":
\co auto mesh = unitSquare();\eco


<a href="#" class="top">top</a>
<hr>
\section Load Load Meshes
\subsection LoadMesh loadMesh
You can use this function to:
\li load a \c .msh file and use the mesh data structure
\li load a \c .geo file and automatically generate a mesh data structure on this geometrical structure

<b>Interface:</b><br>
\co
mesh_ptrtype loadMesh(_mesh, _filename, _refine, _update, physical_are_elementary_regions);
\eco
Required Parameters:
\li \c _mesh: a mesh data structure.

Optional Parameters:
 -  \c _hsize (double): characteristic size of the mesh. This option will edit the \c .geo file and change the variable \c h if defined
   - Default: \c 0.1
   -  Option: \c gmsh.hsize
 -  \c _geo_variables_list (string): Set a list of variable that may be defined in a \c .geo file
   - Default: \c ""
   -  Option: \c gmsh.geo-variables-list
 -  \c _filename (string): filename with extension.
   - Default: \c "feel.geo"
   -  Option: \c gmsh.filename
 -  \c _depends (string): list of files (separated by , or ;) on which \c gmsh.filename depends
   - Default: \c ""
   -  Option: \c gmsh.depends
 - \c _refine (boolean): optionally refine with \p refine levels the mesh.
   - Default: \c 0.
   - Option: \c gmsh.refine
 - \c _update (integer): update the mesh data structure (build internal faces and edges).
   - Default: \c true
 -  \c _physical_are_elementary_regions (boolean): to load specific meshes formats.
   - Default: \c false.
   - Option: gmsh.physical_are_elementary_regions
 - \c _straighten (boolean): in case of curvilinear elements, straighten the elements
   which are not touching with a face the boundary of the domain
   - Default: \c true
   - Option: \c gmsh.straighten
 - \c _partitioner (integer): define the mesh partitioner to use:
   - Default: \c 1 (if Metis is available) \c 0 if not (CHACO)
   - Option: gmsh.partitioner


<br>
The file you want to load has to be in an appropriate repository.<br>
\feel looks for \c .geo and \c .msh files in the following directories (in this order):
\li current path
\li paths that went through changeRepository(), it means that we look for example into the path from which the executable was run
\li localGeoRepository() which is usually \c "$HOME/feel/geo"  (cf: \ref Environment )
\li systemGeoRepository() which is usually \c "$FEELPP_DIR/share/feel/geo" (cf: \ref Environment)


<b>Examples:</b><br>
Load a mesh data structure from the file \c "$HOME/feel/mymesh.msh".
\co
auto mesh = loadMesh(_mesh=new mesh_type,
                     _filename="mymesh.msh");
\eco
<br>
Load a geometric structure from the file \c "./mygeo.geo" and automatically create a mesh data structure.
\co
auto mesh = loadMesh(_mesh=new mesh_type,
                     _filename="mygeo.geo");
\eco
<br>
Create a mesh data structure from the file \c "./feel.geo".
\co
auto mesh = loadMesh(_mesh=new Mesh<Simplex< 2 > > );
\eco

\subsection LoadGMSh loadGMSHMesh
In order to load only \c .msh file, you can also use the loadGMSHMesh.

<b>Interface:</b><br>
\co
mesh_ptrtype loadGMSHMesh(_mesh, _filename, _refine, _update, _physical_are_elementary_regions);
\eco
Required Parameters:
\li \c _mesh: a mesh data structure.
\li \c _filename: filename with extension.

Optional Parameters:
\li \c _refine: optionally refine with \p refine levels the mesh. Default =\c 0.
\li \c _update: update the mesh data structure (build internal faces and edges). Default =\c true.
\li \c _physical_are_elementary_regions: to load specific meshes formats. Default = \c false.

The file you want to load has to be in an appropriate repository. See \ref LoadMesh.

<b>Examples:</b><br>
From \c "doc/manual/heatns.cpp":
\co
 mesh_ptrtype mesh = loadGMSHMesh( _mesh=new mesh_type,
                                   _filename="piece.msh",
                                   _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER );
\eco

From \c "applications/check/check.cpp":
\co
mesh = loadGMSHMesh( _mesh=new mesh_type,
                     _filename=soption("filename"),
                     _rebuild_partitions=(Environment::worldComm().size() > 1),
                     _update=MESH_RENUMBER|MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK );
\eco


\section Create Create Meshes
\subsection CreateGMSHMesh createGMSHMesh
<b>Interface:</b><br>
\co
mesh_ptrtype createGMSHMesh(_mesh, _desc, _h, _order, _parametricnodes, _refine, _update, _force_rebuild, _physical_are_elementary_regions);
\eco
Required Parameters:
\li \c _mesh: mesh data structure.
\li \c _desc: descprition. See further.

Optional Parameters:
\li \c _h: characteristic size. Default = \c 0.1.
\li \c _order: order. Default = \c 1.
\li \c _parametricnodes: Default = \c 0.
\li \c _refine: optionally refine with \p refine levels the mesh. Default =\c 0.
\li \c _update: update the mesh data structure (build internal faces and edges). Default =\c true.
\li \c _force_rebuild: rebuild mesh if already exists. Default = \c false.
\li \c _physical_are_elementary_regions: to load specific meshes formats. Default = \c false.

To generate your mesh you need a description parameter. This one can be create by one the two following function.

\subsection Geo geo
Use this function to create a description from a \c .geo file.

\Interface
\co
gmsh_ptrtype geo(_filename, _h, _dim, _order, _files_path);
\eco

Required Parameters:
\li \c filename: file to load.

Optional Parameters:
\li \c _h: characteristic size of the mesh. Default = \c 0.1.
\li \c _dim: dimension. Default = \c 3.
\li \c _order: order. Default = \c 1.
\li \c _files_path: path to the file. Default = \c localGeoRepository().

The file you want to load has to be in an appropriate repository. See \ref LoadMesh.

\Examples
From \c "doc/manual/heat/ground.cpp":
\co
mesh = createGMSHMesh( _mesh=new mesh_type,
                       _desc=geo( _filename="ground.geo",
                                  _dim=2,
                                  _order=1,
                                  _h=meshSize ) );
\eco

From \c "doc/manual/fd/penalisation.cpp":
\co
mesh = createGMSHMesh( _mesh=new mesh_type,
                       _desc=geo( _filename=File_Mesh,
                                  _dim=Dim,
                                  _h=Environment::vm(_name="hsize").template as<double>() ),
                                  _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER );
\eco


\subsection Domain domain
Use this function to generate a simple geometrical domain from parameters.

\Interface
\co
gmsh_ptrtype domain(_name, _shape, _h, _dim, _order, _convex, \
                    _addmidpoint, _xmin, _xmax, _ymin, _ymax, _zmin, _zmax);
\eco

Required Parameters:
\li \c _name: name of the file that will ge generated without extension.
\li \c _shape: shape of the domain to be generated (simplex or hypercube).

Optional Parameters:
\li \c _h: characteristic size of the mesh. Default = \c 0.1.
\li \c _dim: dimension of the domain. Default = \c 2.
\li \c _order: order of the geometry. Default = \c 1.
\li \c _convex: type of convex used to mesh the domain. Default = \c simplex.

\li \c _addmidpoint: add middle point. Default = \c true.
\li \c _xmin: minimum x coordinate. Default = \c 0.
\li \c _xmax: maximum x coordinate. Default = \c 1.
\li \c _ymin: minimum y coordinate. Default = \c 0.
\li \c _ymax: maximum y coordinate. Default = \c 1.
\li \c _zmin: minimum z coordinate. Default = \c 0.
\li \c _zmax: maximum z coordinate. Default = \c 1.

\Examples
From \c "doc/manual/laplacian/laplacian.ccp":
\co
mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                    _desc=domain( _name=( boost::format( "%1%-%2%" ) % shape % Dim ).str() ,
                                                  _usenames=true,
                                                  _shape=shape,
                                                  _h=meshSize,
                                                  _xmin=-1,
                                                  _ymin=-1 ) );
\eco

From \c "doc/manual/stokes/stokes.cpp":
\co
mesh = createGMSHMesh( _mesh=new mesh_type,
                       _desc=domain( _name=(boost::format("%1%-%2%-%3%")%"hypercube"%convex_type().dimension()%1).str() ,
                                     _shape="hypercube",
                                     _dim=convex_type().dimension(),
                                     _h=meshSize ) );
\eco

From \c "doc/manual/solid/beam.cpp":
\co
mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                    _update=MESH_UPDATE_EDGES|MESH_UPDATE_FACES|MESH_CHECK,
                                    _desc=domain( _name=( boost::format( "beam-%1%" ) % nDim ).str(),
                                                  _shape="hypercube",
                                                  _xmin=0., _xmax=0.351,
                                                  _ymin=0., _ymax=0.02,
                                                  _zmin=0., _zmax=0.02,
                                                  _h=meshSize ) );
\eco


<a href="#" class="top">top</a>
<hr>
\section MeshTodo Todo
\co
straightenMesh
\eco


<a href="#" class="top">top</a>
<hr>
\li \b Next: \ref Spaces
*/
}
