Examples
========

\section example-quickstart Examples: Quickstart

\subsection example-cg-laplacian Hello, Laplacian

\snippet cg_laplacian.hpp cg_laplacian

\subsection example-qs-stokes Hello, Stokes

First we define the environment

\snippet qs_stokes.cpp stokes-env

then the mesh

\snippet qs_stokes.cpp stokes-mesh

then the function solving the Stokes equations

\snippet qs_stokes.cpp stokes

We can now call this functions with different velocity/pressure function spaces: 
- P1P0, P1P0 which aren't inf-sup stable
- P2P1 (aka TaylorHood) which is inf-sup stable

\snippet qs_stokes.cpp stokes-space
