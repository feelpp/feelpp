
/* -*- mode: c++; coding: utf-8 -*- */
namespace Feel {
/**
\page Integrals Integrations, Operators and Norms
\ingroup grIntegrals

\tableofcontents

\li \b Previous: \ref Spaces
\li \b Next: \ref Forms

<hr>
\section Integrals_Introduction Introduction
You should be able to create a mesh now. If it is not the case, get back to the section \ref Mesh.<br>

To use the tools of this sections, you have to precise the domain range using the following keywords:
<table class="manual">
<tr><th>\feel Keyword</th><th>Description</th></tr>
<tr><td> \co elements(mesh) \eco</td><td>All the elements of a mesh</td></tr>
<tr><td> \co markedelements(mesh, id) \eco</td><td>The precise element defined by the id.<br>It can be any element (line, surface, domain, and so on).</td></tr>
<tr><td> \co faces(mesh) \eco</td><td>All the faces of the mesh.</td></tr>
<tr><td> \co markedfaces(mesh) \eco</td><td>All the faces of the mesh which are marked.</td></tr>
<tr><td> \co boundaryfaces(mesh) \eco</td><td>All elements that own a topological dimension one below the mesh. <br>For example, if you mesh is a 2D one, \c boundaryfaces(mesh) will return all the lines (because of dimension \f$2-1=1\f$).<br>These elements which have one dimension less, are corresponding to the boundary faces.</td></tr>
<tr><td> \co internalelements(mesh) \eco</td><td>All the elements of the mesh which are stricly within the domain that is to say they do not share a face with the boundary.</td></tr>
<tr><td> \co boundaryelements(mesh) \eco</td><td>All the elements of the mesh which share a face with the boundary of the mesh.</td></tr>
<tr><td> \co edges(mesh) \eco</td><td>All the edges of the mesh.</td></tr>
<tr><td> \co boundaryedges(mesh) \eco</td><td>All boundary edges of the mesh.</td></tr>
</table>


<a href="#" class="top">top</a>
<hr>
\section Integral Integrals
\subsection Integrate integrate
Thank to its finite element embedded language, \feel has its owned <tt>integrate()</tt> function.

\Interface
\co
  integrate( _range, _expr, _quad, _geomap );
\eco
please notice that the order of the parameter is not important, these are <tt>boost</tt> parameters, so you can enter them in the order you want. <br>
To make it clear, there are two required parameters and 2 optional and they of course can be entered in any order
provided you give the parameter name. If you don't provide the parameter name (that is to say \c _range= or the others) they must be entered in the order they are described
below.

Required parameters:
\li <tt>_range</tt>  = domain of integration
\li <tt>_expr</tt>  = integrand expression

Optional parameters:
\li <tt>_quad</tt>  = quadrature to use instead of the default one, wich means <tt>_Q<integer>()</tt> where the integer is the polynomial order to integrate exactely
\li <tt>_geomap</tt>  = type of geometric mapping to use, that is to say:
<table class="manual">
<tr><th>Feel Parameter</th><th>Description</th></tr>
<tr><td>\co GEOMAP_HO\eco</td><td>High order approximation (same of the mesh)</td></tr>
<tr><td>\co GEOMAP_OPT\eco</td><td>Optimal approximation:<br> high order on boundary elements<br> order 1 in the interior</td></tr>
<tr><td>\co GEOMAP_01\eco</td><td>Order 1 approximation (same of the mesh)</td></tr>
</table>

\Examples
From \c "doc/manual/tutorial/dar.cpp":
\co
  form1( ... ) = integrate( _range = elements( mesh ),
                            _expr = f*id( v ) );
\eco

From \c "doc/manual/tutorial/myintegrals.cpp":
\co
  // compute integral f on boundary
  double intf_3 = integrate( _range = boundaryfaces( mesh ),
                             _expr = f );
\eco

From \c "doc/manual/advection/advection.cpp":
\co
  form2( _test = Xh, _trial = Xh, _matrix = D ) +=
    integrate( _range = internalfaces( mesh ),
               _quad = _Q<2*Order>(),
               _expr = ( averaget( trans( beta )*idt( u ) ) * jump( id( v ) ) )
                     + penalisation*beta_abs*( trans( jumpt( trans( idt( u ) ) ) )*jump( trans( id( v ) ) ) ),
               _geomap = geomap );
\eco

From \c "doc/manual/laplacian/laplacian.cpp":
\co
 auto l = form1( _test=Xh, _vector=F );
 l = integrate( _range = elements( mesh ),
                _expr=f*id( v ) ) +
     integrate( _range = markedfaces( mesh, "Neumann" ),
                _expr = nu*gradg*vf::N()*id( v ) );
\eco


\subsection Integrals_Computing Computing my first Integrals
This part explains how to integrate on a mesh with \feel (source \c "doc/manual/tutorial/myintegrals.cpp").

Let's consider the domain \f$\Omega=[0,1]^d\f$ and associated meshes.<br>
Here, we want to integrate the following function
<br><center>\f$
\begin{aligned}
f(x,y,z) = x^2 + y^2 + z^2
\end{aligned}
\f$</center><br>
on the whole domain \f$\Omega\f$ and on part of the boundary \f$\Omega\f$.

There is the appropriate code:
\co
int
main( int argc, char** argv )
{
    // Initialize Feel++ Environment
    Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _about=about( _name="myintegrals" ,
                                   _author="Feel++ Consortium",
                                   _email="feelpp-devel@feelpp.org" ) );

    // create the mesh (specify the dimension of geometric entity)
    auto mesh = unitHypercube<3>();

    // our function to integrate
    auto f = Px()*Px() + Py()*Py() + Pz()*Pz();

    // compute integral of f (global contribution)
    double intf_1 = integrate( _range = elements( mesh ),
                               _expr = f ).evaluate()( 0,0 );

    // compute integral of f (local contribution)
    double intf_2 = integrate( _range = elements( mesh ),
                               _expr = f ).evaluate(false)( 0,0 );

    // compute integral f on boundary
    double intf_3 = integrate( _range = boundaryfaces( mesh ),
                               _expr = f ).evaluate()( 0,0 );

    std::cout << "int global ; local ; boundary" << std::endl
              << intf_1 << ";" << intf_2 << ";" << intf_3 << std::endl;
}
\eco


<a href="#" class="top">top</a>
<hr>
\section Operators Operators
\subsection Project project
It is also possible to make projections with the library.

\Interface
\co
  project( _range, _space, _expr, _geomap );
\eco

Required parameters:
\li <tt>_space</tt>: the space in which lives the projected expression, it should be a nodal function space
\li <tt>_expr</tt>: the expression to project

Optional parameters:
\li <tt>_range</tt>: the domain for the projection. Default = all elements from <tt>space->mesh()</tt>
\li <tt>_geomap</tt>: type of geometric mapping. Default = <tt>GEOMAP_OPT</tt>

\Examples
From \c "doc/manual/laplacian/laplacian.cpp":
\co
  element_type e( Xh, "e" );
  e = project( _space = Xh,
               _range = elements( mesh ),
               _expr = g );
\eco

From \c "doc/manual/heatns/convection_run.cpp":
\co
tn = project( _space = Xh->functionSpace<2>(),
              _range = elements( mesh ),
              _expr = constant( 300 ) );
\eco



\subsection Mean mean
Let \f$f\f$ a bounded function on domain \f$\Omega\f$. You can evaluate the mean value:
<br><center>\f$
  \begin{aligned}
 \bar{f}&=\frac{1}{|\Omega|}\int_\Omega f\\
&=\frac{1}{\int\limits_\Omega 1}\int_\Omega f.
  \end{aligned}
\f$</center><br>

\Interface
\co
  mean( _range, _expr, _quad, _geomap );
\eco

Required parameters:
\li <tt>_range</tt> = domain of integration
\li <tt>_expr</tt> = mesurable function

Optional parameters:
\li <tt>_quad</tt> = quadrature to use. Default = \lstinline!_Q<integer>()!
\li <tt>_geomap</tt> = type of geometric mapping. Default = <tt>GEOMAP_OPT</tt>

\Examples
From \c "doc/manual/stokes/stokes.cpp":
\snippet stokes.cpp mean

<a href="#" class="top">top</a>
<hr>
\section Norms Norms
\subsection NormL2 normL2
Let \f$f \in L^2(\Omega)\f$ you can evaluate the L2 norm:
<br><center>\f$
  \begin{aligned}
\parallel f\parallel_{L^2(\Omega)}=\sqrt{\int_\Omega |f|^2}
  \end{aligned}
\f$</center><br>

\Interface
\co
  normL2( _range, _expr, _quad, _geomap );
\eco
or squared norm:
\co
  normL2Squared( _range, _expr, _quad, _geomap );
\eco

Required parameters:
\li <tt>_range</tt> = domain of integration
\li <tt>_expr</tt>  = mesurable function

Optional parameters:
\li <tt>_quad</tt>  = quadrature to use. Default = <tt>_Q<integer>()</tt>
\li <tt>_geomap</tt>  = type of geometric mapping. Default = <tt>GEOMAP_OPT</tt>

\Examples
From \c "doc/manual/laplacian/laplacian.cpp":
\co
  double L2error =normL2( _range=elements( mesh ),
                          _expr=( idv( u )-g ) );
\eco

From \c "doc/manual/stokes/stokes.cpp":
\snippet stokes.cpp norml2

\subsection NormH1 normH1
In the same idea, you can evaluate the H1 norm or semi norm, for any function \f$f \in H^1(\Omega)\f$:
<br><center>\f$
\begin{aligned}
 \parallel f \parallel_{H^1(\Omega)}&=\sqrt{\int_\Omega |f|^2+|\nabla f|^2}\\
&=\sqrt{\int_\Omega |f|^2+\nabla f*\nabla f^T}\\
|f|_{H^1(\Omega)}&=\sqrt{\int_\Omega |\nabla f|^2}
\end{aligned}
\f$</center><br>

\Interface
\co
  normH1( _range, _expr, _grad_expr, _quad, _geomap );
\eco
or semi norm:
\co
  normSemiH1( _range, _grad_expr, _quad, _geomap );
\eco

Required parameters:
\li <tt>_range</tt> = domain of integration
\li <tt>_expr</tt> = mesurable function
\li <tt>_grad_expr</tt> = gradient of function (Row vector!)

Optional parameters:
\li <tt>_quad</tt> = quadrature to use. Default = <tt>_Q<integer>()</tt>
\li <tt>_geomap</tt> = type of geometric mapping. Default = <tt>GEOMAP_OPT</tt>


\Examples
With expression:
\co
  auto g = sin(2*pi*Px())*cos(2*pi*Py());
  auto gradg = 2*pi*cos(2* pi*Px())*cos(2*pi*Py())*oneX() \
	           -2*pi*sin(2*pi*Px())*sin(2*pi*Py())*oneY();
// There gradg is a column vector!
// Use trans() to get a row vector
  double normH1_g = normH1( _range=elements(mesh),\
  			                _expr=g,\
			                _grad_expr=trans(gradg) );
\eco
With test or trial function \c u:
\co
  double errorH1 = normH1( _range=elements(mesh),\
  			               _expr=(u-g),\
    			           _grad_expr=(gradv(u)-trans(gradg)) );
\eco



\subsection normLinf normLinf
Let \f$f\f$ a bounded function on domain \f$\Omega\f$. You can evaluate the infinity norm:
<br><center>\f$
\parallel f \parallel_\infty=\sup_\Omega(|f|)
\f$</center><br>

\Interface:
\co
  normLinf( _range, _expr, _pset, _geomap );
\eco

Required parameters:
\li <tt>_range</tt> = domain of integration
\li <tt>_expr</tt> = mesurable function
\li <tt>_pset</tt> = set of points (e.g. quadrature points)

Optional parameters:
\li <tt>_geomap</tt> = type of geometric mapping. Default = <tt>GEOMAP_OPT</tt>


\Examples
\co
  auto uMax = normLinf( _range=elements(mesh),\
		        _expr=idv(u),\
		        _pset=_Q<5>() );
\eco



<a href="#" class="top">top</a>
<hr>

\li \b Next: \ref Solver
*/
}
