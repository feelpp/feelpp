/* -*- mode: c++; coding: utf-8 -*- */
namespace Feel {
/*! \page TutorialIntegrals Computing integrals over mesh

\tableofcontents

The next step is to compute integrals over the mesh. The source code is available in \c myintegrals.cpp.

\section TutorialIntegralsSteps Step by step explanations

We start by loading a Mesh in 2D
\snippet myintegrals.cpp mesh

then we define the Feel++ expression that we are going to integrate using the \c
soption function that retrieves the command line option string \c functions.g. We then transform this string into a Feel++ expression using \c expr().

\snippet myintegrals.cpp expression

then We compute two integrals over the domain and its boundary respectively

 @li \f$\int_\Omega g\f$

 @li \f$\int_{\partial \Omega} g\f$

and we print the results to the screen.

\snippet myintegrals.cpp integrals

\remark Only the rank 0 process (thanks to \c Environment::isMasterRank() prints
to the screen as the result is the same over all mpi processes if the
application was run in parallel. Note also that the code actually prints the
expression passed by the user through the command line option \c functions.g.

\section TutorialIntegralsResults Some results

We start with the following function \f$g=1\f$. Recall that by default the
domain is the unit square, hence the \f$\int_\Omega g\f$ and \f$\int_{\partial
\Omega} g\f$ should be equal to 1 and 4 respectively.

\code{.sh}
shell> ./feelpp_doc_myintegrals --functions.g=1
int_Omega 1 = 1
int_{boundary of Omega} 1 = 4
\endcode

Now we try \f$g=x\f$. We need to tell Feel++ what are the symbols associated
with the expression: here the symbol is \c x and it works as follows

\code{.sh}
shell> ./feelpp_doc_myintegrals --functions.g=x:x
int_Omega x = 0.5
int_{boundary of Omega} x = 2
\endcode

Recall that \c : is a separator between the expression and each symbol composing it.

Now we try \f$g=x y\f$. We need to tell Feel++ what are the symbols associated
with the expression: here the symbol is \c x and \c y and it works as follows

\code{.sh}
shell> ./feelpp_doc_myintegrals --functions.g=x*y:x:y
int_Omega x*y = 0.25
int_{boundary of Omega} x*y = 1
\endcode

More complicated functions are of course doable, such as \f$g=\sin( x y )\f$.

\code{.sh}
./feelpp_doc_myintegrals --functions.g="sin(x*y):x:y"
int_Omega sin(x*y) = 0.239812
int_{boundary of Omega} sin(x*y) = 0.919395
\endcode

Here is the last example in parallel over 4 processors which returns, of
course, the exact same results as in sequential

\code{.sh}
mpirun -np 4 ./feelpp_doc_myintegrals --functions.g="sin(x*y):x:y"
int_Omega sin(x*y) = 0.239812
int_{boundary of Omega} sin(x*y) = 0.919395
\endcode

Finally we can change the type of domain and compute the area and perimeter of the unit disk as follows
\code{.sh}
./feelpp_doc_myintegrals --functions.g="1:x:y" --gmsh.domain.shape=ellipsoid --gmsh.hsize=0.05
int_Omega 1 = 0.784137
int_{boundary of Omega} 1 = 3.14033
\endcode

\remark Note that we don't get the exact results due to the fact that
\f$\Omega_h = \cup_{K \in \mathcal{T}_h} K\f$ which we use for the numerical
integration is different from the exact domain \f$\Omega = \{ (x,y)\in
\mathbb{R}^2 | x^2+y^2 < 1\}\f$.

\section TutorialIntegralsCode Complete code

The complete code reads as follows

\snippet myintegrals.cpp all

to compile just type
\verbatim
make feelpp_doc_myintegrals
\endverbatim

*/
}
