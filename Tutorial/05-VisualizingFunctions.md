/* -*- mode: c++; coding: utf-8 -*- */
namespace Feel {
/*! \page TutorialVisualize Visualizing functions over a mesh


The next step is to visualize function over the mesh. The source code is
available in \c myexporter.cpp.

@par Loading a Mesh in 2D
Here, we generate a second order mesh,
\snippet myexporter.cpp mesh
and one of first order.
\snippet myexporter.cpp P1_mesh


@par Constructing a function space
here, we generate a second order function space,
\snippet myexporter.cpp space

@par Defining a (scalar) function over the function space
\snippet myexporter.cpp function

@par Exporter
We create now three exporter:
- once generated on the P1 mesh
- once generated on the P2 mesh
- once generated on a P1 mesh extracted from P2 mesh

\snippet myexporter.cpp exporter

@par Adding function to save
Here we save the function many times.
That is here not relevant but you may want to simulate process over time. 
\snippet myexporter.cpp adding

@par Actually saving 
\snippet myexporter.cpp save


\section TutorialExprCode Complete code

The complete code reads as follows

\snippet myexporter.cpp all

to compile just type
\verbatim
make feelpp_doc_myexporter
\endverbatim
to execute just type
\verbatim
./feelpp_doc_myexporter
\endverbatim


@par Reading saved data
You can visualize data via
- [ensight](https://www.ceisoftware.com/),
- [paraview](www.paraview.org/),
- [gmsh](http://geuz.org/gmsh).
 
The results files are in \c $HOME/feel/myexporter/np_1 or \c $FEELPP_WORKDIR/feel/myexporter/np_1

*/
}
