## Introduction

[Feel++](www.feelpp.org) is a C++ library for arbitrary order Galerkin methods (e.g. finite and spectral element methods) continuous or discontinuous in 1D 2D and 3D. The objectives of this framework is quite ambitious; ambitions which could be express in various ways such as :

  - the creation of a versatile mathematical kernel solving easily problems using different techniques thus allowing testing and comparing methods, e.g. cG versus dG,
  - the creation of a small and manageable library which shall nevertheless encompass a wide range of numerical methods and techniques,
  - build mathematical software that follows closely the mathematical abstractions associated with partial differential equations (PDE)
  - the creation of a library entirely in C++ allowing to create C++ complex and typically multi-physics applications such as fluid-structure interaction or mass transport in haemodynamic


Some basic installation procedure are available in the [INSTALL](INSTALL.md) file, the detailled process is available [here](http://www.feelpp.org/docs/develop/BuildingP.html)

## Releases

Here are the latest releases of Feel++ 

   - Feel++ [0.99.0](https://github.com/feelpp/feelpp/releases/tag/v0.99.0-final.1) [![DOI](https://zenodo.org/badge/3881/feelpp/feelpp.png)](http://dx.doi.org/10.5281/zenodo.11624)

## Build Status

Feel++ uses Travis-CI for continuous integration.
Travis-CI Build Status :

  - develop branch : [![Build Status](https://travis-ci.org/feelpp/feelpp.svg?branch=develop)](https://travis-ci.org/feelpp/feelpp)
  - master branch : [![Build Status](https://travis-ci.org/feelpp/feelpp.svg?branch=master)](https://travis-ci.org/feelpp/feelpp)

## Documentation

  - develop branch : [Feel++ Online Reference Manual](http://www.feelpp.org/docs/develop/index.html)
  - master branch (latest release) : [Feel++ Online Reference Manual](http://www.feelpp.org/docs/master/index.html)

## Features

  - 1D 2D and 3D (including high order) geometries and also lower topological dimension 1D(curve) in 2D and 3D or 2D(surface) in 3D
  - continuous and discontinuous arbitrary order Galerkin Methods in 1D, 2D and 3D including finite and spectral element methods
  - domain specific embedded language in C++ for variational formulations
  - interfaced with [PETSc](http://www.mcs.anl.gov/petsc/) for linear and non-linear solvers
  - seamless parallel computations using PETSc
  - interfaced with [SLEPc](http://www.grycap.upv.es/slepc/) for large-scale sparse standard and generalized eigenvalue  solvers
  - supports [Gmsh](http://www.geuz.org/gmsh) for mesh generation
  - supports [Gmsh](http://www.geuz.org/gmsh) for post-processing (including on high order geometries)
  - supports [Paraview](http://www.paraview.org) for post-processing


## Examples

### Laplacian in 2D using P3 Lagrange basis functions

Here is a full example to solve
$$-\Delta u = f \mbox{ in } \Omega,\quad u=g \mbox{ on } \partial \Omega$$

```cpp
#include <feel/feel.hpp>

int main(int argc, char**argv )
{
    using namespace Feel;
	Environment env( _argc=argc, _argv=argv,
                     _desc=feel_options(),
                     _about=about(_name="qs_laplacian",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    auto mesh = unitSquare();
    auto Vh = Pch<1>( mesh );
    auto u = Vh->element();
    auto v = Vh->element();

    auto l = form1( _test=Vh );
    l = integrate(_range=elements(mesh),
                  _expr=id(v));

    auto a = form2( _trial=Vh, _test=Vh );
    a = integrate(_range=elements(mesh),
                  _expr=gradt(u)*trans(grad(v)) );
    a+=on(_range=boundaryfaces(mesh), _rhs=l, _element=u,
          _expr=constant(0.) );
    a.solve(_rhs=l,_solution=u);

    auto e = exporter( _mesh=mesh, _name="qs_laplacian" );
    e->add( "u", u );
    e->save();
    return 0;
}
```


### Bratu equation in 2D

Here is a full non-linear example - the Bratu equation - to solve
$$-\Delta u + e^u = 0 \mbox{ in } \Omega,\quad u=0 \mbox{ on } \partial \Omega$$.

```cpp
#include <feel/feel.hpp>

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description bratuoptions( "Bratu problem options" );
    bratuoptions.add_options()
    ( "lambda", Feel::po::value<double>()->default_value( 1 ),
                "exp() coefficient value for the Bratu problem" )
    ( "penalbc", Feel::po::value<double>()->default_value( 30 ),
                 "penalisation parameter for the weak boundary conditions" )
    ( "hsize", Feel::po::value<double>()->default_value( 0.1 ),
               "first h value to start convergence" )
    ( "export-matlab", "export matrix and vectors in matlab" )
    ;
    return bratuoptions.add( Feel::feel_options() );
}

/**
 * Bratu Problem
 *
 * solve \f$ -\Delta u + \lambda \exp(u) = 0, \quad u_\Gamma = 0\f$ on \f$\Omega\f$
 */
int
main( int argc, char** argv )
{

    using namespace Feel;
	Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=about(_name="bratu",
                                  _author="Christophe Prud'homme",
                                  _email="christophe.prudhomme@feelpp.org"));
    auto mesh = unitSquare();
    auto Vh = Pch<3>( mesh );
    auto u = Vh->element();
    auto v = Vh->element();
    double penalbc = option(_name="penalbc").as<double>();
    double lambda = option(_name="lambda").as<double>();

    auto Jacobian = [=](const vector_ptrtype& X, sparse_matrix_ptrtype& J)
        {
            auto a = form2( _test=Vh, _trial=Vh, _matrix=J );
            a = integrate( elements( mesh ), gradt( u )*trans( grad( v ) ) );
            a += integrate( elements( mesh ), lambda*( exp( idv( u ) ) )*idt( u )*id( v ) );
            a += integrate( boundaryfaces( mesh ),
               ( - trans( id( v ) )*( gradt( u )*N() ) - trans( idt( u ) )*( grad( v )*N()  + penalbc*trans( idt( u ) )*id( v )/hFace() ) );
        };
    auto Residual = [=](const vector_ptrtype& X, vector_ptrtype& R)
        {
            auto u = Vh->element();
            u = *X;
            auto r = form1( _test=Vh, _vector=R );
            r = integrate( elements( mesh ), gradv( u )*trans( grad( v ) ) );
            r +=  integrate( elements( mesh ),  lambda*exp( idv( u ) )*id( v ) );
            r +=  integrate( boundaryfaces( mesh ),
               ( - trans( id( v ) )*( gradv( u )*N() ) - trans( idv( u ) )*( grad( v )*N() ) + penalbc*trans( idv( u ) )*id( v )/hFace() ) );
        };
    u.zero();
    backend()->nlSolver()->residual = Residual;
    backend()->nlSolver()->jacobian = Jacobian;
    backend()->nlSolve( _solution=u );

    auto e = exporter( _mesh=mesh );
    e->add( "u", u );
    e->save();
}
```
