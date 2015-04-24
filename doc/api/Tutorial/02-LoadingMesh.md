Loading a Mesh {#TutorialMesh}
==================

\tableofcontents


The next step is to load a mesh. The source code is available in \c mymesh.cpp.

# TutorialMeshSteps Step by step {#explanations}

## Loading a Mesh in 2D {#first}
The `loadMesh` function has a `_name` option set by default as the default value of the `--gmsh.filename` option that point either to a `.geo`, either to a `.msh` file.

\snippet mymesh.cpp load

## Exporting the Mesh for visualisation {#visu}

Please refere [here](TutorialVisualize.html) for explanations.
\snippet mymesh.cpp export

# Complete code {#TutorialMeshCode}

\snippet mymesh.cpp all

