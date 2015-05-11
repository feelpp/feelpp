/* -*- mode: c++; coding: utf-8 -*- */
namespace Feel {
/*! \page Forms Forms and Solver


\tableofcontents


<hr>

We suppose in this section that you know how to define your mesh (see \ref Mesh) and your function spaces (see \ref Spaces). You may need integration tools too (see \ref Integrals).

There are \feel tools you need to create linear and bilinear forms in order to solve variational formulation.

Notations:
\li \c u: element from your trial function space (unknown function)
\li \c v: element from your test function space

\section Forms_Building Building Forms
\subsection form1 form1
\Interface
\co
form1(_test, _init);
\eco
Required Parameters:
\li \c _test: test function space.

Optional Parameters:
\li \c _init: Default = \c false.

By default, a new linear form is:
\f[
l(v)=\int_\Omega v
\f]
Then you can customize it using integration tools.

\Examples
From `mylaplacian.cpp`
\snippet mylaplacian.cpp marker_form1

From `myadvection.cpp`
\snippet myadvection.cpp marker_form1

Notice that \c += operator is working with linear and bilinear forms.


\subsection form2 form2
\Interface
\co
form2(_trial, _test, _init);
\eco
Required Parameters:
\li \c _trial: test function space
\li \c _test: trial function space

Optional Parameters:
\li \c _init: Default = \c false.

By default, a new bilinear form is:
\f[
a(u,v)=\int_\Omega uv
\f]
Then you can custom it using integrations tools

\Example
From `mylaplacian.cpp`
\snippet mylaplacian.cpp marker_form2

From `mystokes.cpp`:
\snippet mystokes.cpp marker_form2

Notice that \c += operator is working with linear and bilinear forms.


<a href="#" class="top">top</a>
<hr>
\section Solver Solver
In this section we present syntax to solve variational formulations. For more general linear problems see \ref Linear.<br>

\subsection solve solve
Once you created your linear and bilinear forms you can use the \c solve() function on your bilinear form.<br>
The \c solve() function presented there is a method from the class \c BilinearForm.<br>
\Interface
\co
solve(_solution, _rhs, _rebuild, _name);
\eco
Required Parameters:
\li \c _solution: the solution.
\li \c _rhs: right hand side. The linear form.

Optional Parameters:
\li \c _rebuild: rebuild the solver matrix. Default = \c false.
\li \c _name: Default = "".

\Examples
From `laplacian.cpp`:
\snippet mylaplacian.cpp marker_solve

\subsection on on
The function \c on() allows you to add conditions to your bilinear form before using the \c solve function.<br>
\Interface
\co
on(_range, _rhs, _element, _expr);
\eco
Required Parameters:
\li \c _range: domain concerned by this condition (see \ref Integrals ).
\li \c _rhs: right hand side. The linear form.
\li \c _element: element concerned.
\li \c _expr: the condition.

This function is used with += operator.

\Examples
From `mylaplacian.cpp`:
\snippet mylaplacian.cpp marker_on
There we add the condition:\f$ u  =  0  \text{ on }\;\partial\Omega \;\f$.

From `mystokes.cpp`:
\snippet mystokes.cpp marker_on

You can also apply boundary conditions using :
 \co
  a+=on(_range=markedfaces(mesh,"top"),_element=u[Component::Y],_rhs=l,_expr=cst(0.))
\eco

<a href="#" class="top">top</a>
<hr>

*/
}
