Visualizing functions over a mesh {#TutorialVisualize}
======================================



The next step is to visualize function over the mesh. The source code is
available in `myexporter.cpp.`

# Loading a Mesh in 2D {#load}

Here, we generate a second order mesh,
!CODEFILE "code/myexporter.cpp" mesh
and one of first order.
!CODEFILE "code/myexporter.cpp" P1_mesh


# Constructing a function space {#fs}

here, we generate a second order function space,
!CODEFILE "code/myexporter.cpp" space

# Defining a (scalar) function over the function space {#scal}

!CODEFILE "code/myexporter.cpp" function

# Exporter {#exp}

We create now three exporter:
- once generated on the P1 mesh
- once generated on the P2 mesh
- once generated on a P1 mesh extracted from P2 mesh

!CODEFILE "code/myexporter.cpp" exporter

# Adding function to save {#add}

Here we save the function many times.
That is here not relevant but you may want to simulate process over time.
!CODEFILE "code/myexporter.cpp" adding

# Actually saving {#save}

!CODEFILE "code/myexporter.cpp" save


#  Complete code {#TutorialExprCode}

The complete code reads as follows

!CODEFILE "code/myexporter.cpp" all

to compile just type
```cpp
make feelpp_doc_myexporter
```
to execute just type
```cpp
./feelpp_doc_myexporter
```


# Reading saved data {#read}

You can visualize data via
- [ensight](https://www.ceisoftware.com/),
- [paraview](www.paraview.org/),
- [gmsh](http://geuz.org/gmsh).

The results files are in \c $HOME/feel/myexporter/np_1 or \c $FEELPP_WORKDIR/feel/myexporter/np_1
