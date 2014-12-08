/* -*- mode: c++; coding: utf-8 -*- */
namespace Feel {
/*! \page Projection Projection



\li \b Previous: \ref Spaces
\li \b Next: \ref Forms

<hr>

\tableofcontents


\section NodalProjection Nodal Projection

An important keyword is \c project which allows to compute the
interpolant of an expression with respect to a nodal function space over a part
of the mesh  or the whole mesh. The interface is as follows
\code
project( _space=<nodal function space in which the interpolant lives>
         _range=<domain range iterators>,
         _expr=<expression to be interpolated>, .... )
\endcode

Here are some examples:

\code
typedef FunctionSpace<Mesh<Simplex<$d$>,
                      bases<Lagrange<1,Vectorial> > > Xhv_type;
auto Xhv = Xhv_type::New( mesh );
// build a piecewise $\P_1$ vectorial function in Xhv containing the
// coordinates of the vertices the mesh.
auto coord = project( _space=Xhv, _range=elements(mesh), _expr=P() );
// compute the x derivative of the coord function
auto dx_coord = project( _space=Xhv, _range=elements(mesh), _expr=dxv(coord) );
auto dy_coord = project( _space=Xhv, _range=elements(mesh), _expr=dyv(coord) );
\endcode

\section ProjectionProjector  Projection Operator

A projection operator is available that supports:

 - \f$L_2\f$ projection
 - \f$H_1\f$ projection

\section ProjectionExpression  Project Specific Expressions

You may want to customise the projected expression, that is to define your own expression.
We have based that mechanism over the GiNaC librarie or the Functor.

The idea with GiNaC is to provide a string which will be parsed to generate a function:

\snippet myexpression.cpp expr

To use the Functor, you have to define in the Feel namespace a struct with some variable.
That struct will provide an \c operator() interface with specific signature.

\snippet myfunctor.cpp functor

and then, you are allowed to define a functor based expression on that way:

\snippet myfunctor.cpp functors_and_proj

<a href="#" class="top">top</a>
<hr>
*/
}
