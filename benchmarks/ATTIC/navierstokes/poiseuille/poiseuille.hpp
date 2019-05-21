/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2011-03-15

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file stokes_poiseuille.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2011-03-15
 */
#include <feel/feelcore/feel.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>

template<int Dim>
void
poiseuille( int argc, char** argv )
{
    double penalbc_u=100;
    double penalbc_p=10;
    double mu=1;
    // Declare the supported options.
    namespace po = boost::program_options;
    po::options_description desc( "Allowed options" );
    desc.add_options()
    ( "help", "produce help message" )
    ( "mu", po::value<double>( &mu )->default_value( 1 ), "viscosity" )
    ( "penalbc_u", po::value<double>( &penalbc_u )->default_value( 100 ), "penalisation u" )
    ( "penalbc_p", po::value<double>( &penalbc_p )->default_value( 10 ), "penalisation p" )
    ;

    po::variables_map vm;
    po::store( po::parse_command_line( argc, argv, desc ), vm );
    po::notify( vm );

    using namespace Feel;
    using namespace Feel::vf;
    Feel::Environment env( argc, argv );
    typedef Mesh<Simplex<2> > mesh_type;
    typedef FunctionSpace<mesh_type, bases<Lagrange<2,Vectorial>, Lagrange<1,Scalar> > > fs_type;

    auto mesh = createGMSHMesh( _mesh=new mesh_type,
                                _desc=domain( _name=( boost::format( "hypercube-%1%" ) % Dim ).str(),
                                        _usenames=false,
                                        _shape="hypercube",
                                        _dim=Dim,
                                        _xmin=0,_xmax=4,
                                        _ymin=0,_ymax=1,
                                        _h=0.1 ) );
    auto Xh = fs_type::New( mesh );
    auto U = Xh->element(); // velocity + pressure function
    auto u = U.element<0>(); // extract velocity
    auto p = U.element<1>(); // extract pressure
    auto V = Xh->element();
    auto v = V.element<0>();
    auto q = V.element<1>();
    std::cout << "mesh and space created" << std::endl;

    //auto deft=.5*(gradt(u)+trans(gradt(u)));
    //auto def=.5*(grad(u)+trans(grad(u)));
    auto deft=gradt( u );
    auto def=grad( u );
    auto SigmaNt = -idt( p )*N() + 2*deft*N();
    auto SigmaN = -id( p )*N() + 2*def*N();


    // Stokes matrix
    auto S = Backend<double>::build(soption("backend"))->newMatrix( Xh, Xh );
    form2( _test=Xh, _trial=Xh, _matrix=S, _init=true ) = integrate( elements( mesh ), trace( deft*trans( def ) ) ) ;
    form2( _test=Xh, _trial=Xh, _matrix=S ) += integrate( elements( mesh ), - idt( p )*div( v ) + id( q )*divt( u ) );
    //form2(_test=Xh, _trial=Xh, _matrix=S ) += integrate( elements(mesh), 1e-6*idt(p)*id(q) );
    using namespace std;
    using namespace boost::assign; // bring 'operator+=()' into scope

    vector<int> markers;
    markers += 1,2,4; // insert values at the end of the container
    BOOST_FOREACH( int marker, markers )
    {
        form2( _test=Xh, _trial=Xh, _matrix=S ) += integrate( markedfaces( mesh,marker ), -trans( SigmaNt )*id( v ) );
        form2( _test=Xh, _trial=Xh, _matrix=S ) += integrate( markedfaces( mesh,marker ), -trans( SigmaN )*idt( v ) );
        form2( _test=Xh, _trial=Xh, _matrix=S ) += integrate( markedfaces( mesh,marker ), penalbc_u*trans( idt( u ) )*id( v )/hFace() );

    }
    std::cout << "bilinear form created" << std::endl;
    // right hand side
    auto F = Backend<double>::build(soption("backend"))->newVector( Xh );

    if ( Dim == 2 )
        form1( _test=Xh, _vector=F,_init=true ) =
            integrate( markedfaces( mesh,1 ),
                       trans( vec( Py()*( 1.-Py() ),cst( 0. ) ) )*( -SigmaN+penalbc_u*id( v )/hFace() ) );

    if ( Dim == 3 )
        form1( _test=Xh, _vector=F,_init=true ) =
            integrate( markedfaces( mesh,1 ),
                       trans( vec( Py()*Pz()*( 1.-Py() )*( 1.-Pz() ),cst( 0. ) ) )*( -SigmaN+penalbc_u*id( v )/hFace() ) );

    std::cout << "linear form created" << std::endl;
    Backend<double>::build(soption("backend"))->solve( _matrix=S, _rhs=F, _solution=U );
    std::cout << "linear system solved" << std::endl;


    auto exporter = Exporter<mesh_type,1>::New( "ensight", ( boost::format( "poiseuille-%1%" ) % Dim ).str() );
    exporter->step( 0 )->setMesh( mesh );
    exporter->step( 0 )->add( "u", u );
    exporter->step( 0 )->add( "p", p );
    exporter->save();
    std::cout << "done" << std::endl;
}

