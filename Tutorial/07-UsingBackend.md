/* -*- mode: c++; coding: utf-8 -*- */
namespace Feel {
/*! \page TutorialBackend Using a backend

\tableofcontents

\section section1 Introduction

After the discretization process, one may have to solve a (non) linear system. Feel++ interfaces with PETSc/SLEPc and Eigen3.
Consider this system
\f[A x = b \f]

We call \c Backend an object that manages the solution strategy to solve it. Some explanation are available \ref Solver and \ref Preconditioner.


Feel++ provides a default backend that is mostly hidden to the final user.
In many examples, you do not have to take care of the backend. You change the backend behavior via the command line or config files.
For example
\verbatim
./feelpp_doc_mybackend --backend.pc-type=id
\endverbatim
will use the identity matrix as a right preconditionner for the default backend.
The size of the preconditionner will be defined from the size of the A matrix.

If you try to solve a different system \f$A_1 y= c\f$ (in size) with the same backend or the default without rebuilding it, it will fail.
\verbatim
backend(_rebuild=true)->solve(_matrix=A1,_rhs=c,_sol=y);
\endverbatim

\remark Each of that options can be retrieved via the \c --help-lib argument in the command line.

\subsection feel3 Non default backend

You may need to manage more than one backend in an application: you have different systems to solve and you want to
keep some already computed objects such as preconditioners.

The default backend is in fact an unnamed backend: in order to distinguish between backend you have to name them. for example
\snippet mybackend.cpp marker_opt

After that, you create the backend object:
\snippet mybackend.cpp marker_obj

\remark Be careful, the backend' name has to match the name you gave at the options step.

Then, you load meshes, creates spaces etc. At solve time, or you solve with the default backend:
\snippet mybackend.cpp marker_default

One of the important backend option is to be able to monitor the residuals and iteration count
\verbatim
./feelpp_doc_mybackend --pc-type=id --ksp-monitor=true --myBackend.ksp-monitor=true
\endverbatim

Finally you can create a named backend:
\snippet mybackend.cpp marker_hm


The whole code example :
\snippet mybackend.cpp marker_main
*/
}
