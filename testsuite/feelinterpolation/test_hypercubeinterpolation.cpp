/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Abdoulaye Samake <abdoulaye.samake1@e.ujf-grenoble.fr>
       Date: 2011-08-20

  Copyright (C) 2008-2010 Universite Joseph Fourier (Grenoble I)

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
   \file test_trace.cpp
   \author Abdoulaye Samake <abdoulaye.samake.@imag.fr>
   \date 2011-08-20
 */

#include <feel/options.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/operatortrace.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>

using namespace Feel;
using namespace Feel::vf;

inline
po::options_description
makeOptions()
{
    po::options_description testoptions( "Test options" );
    testoptions.add_options()
        ( "hsize", po::value<double>()->default_value( 0.1 ), "mesh size" )
        ;
    return testoptions.add( Feel::feel_options() );
}


template<int Dim, int Order>
class Test
    :
public Simget
{
    typedef Simget super;
public:

    typedef Hypercube<Dim,1,Dim> convex_type;
    typedef Mesh< convex_type > mesh_type;
    typedef FunctionSpace<mesh_type, bases<Lagrange<Order,Scalar> > > space_type;

    /**
     * Constructor
     */

    Test()
        :
        super(),
        backend( backend_type::build( this->vm() ) ),
        meshSize( this->vm()["hsize"].template as<double>() )
    {}

    void run();

private:

    backend_ptrtype backend;
    double meshSize;

}; // Test

template<int Dim, int Order>
void
Test<Dim,Order>::run()
{

    std::cout << "------------------------------------------------------------\n";
    std::cout << "Execute Test<" << Dim << "," << Order << ">\n";

    if ( !this->vm().count( "nochdir" ) )
        Environment::changeRepository( boost::format( "testsuite/feeldiscr/%1%/%2%-%3%/P%4%G%4%/h_%5%/" )
                                       % this->about().appName()
                                       % convex_type::name()
                                       % Dim
                                       % Order
                                       % meshSize );

    auto mesh = createGMSHMesh( _mesh=new mesh_type,
                                _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                                _desc=domain( _name=( boost::format( "hypercube-%1%" ) % Dim ).str(),
                                              _addmidpoint=false,
                                              _usenames=false,
                                              _shape="hypercube",
                                              _dim=Dim,
                                              _h=meshSize,
                                              _convex="Hypercube",
                                              //_convex=convex_type::name(),
                                              _xmin=0.,
                                              _xmax=1.,
                                              _ymin=0.,
                                              _ymax=1.,
                                              _substructuring=true
                                              ),
                                _structured=2);


    auto Xh = space_type::New(_mesh=mesh);
    auto TXh = Xh->trace( markedfaces( mesh,"NORTH" ) ) ;

    auto opI = opInterpolation( _domainSpace=Xh,
                                _imageSpace=TXh,
                                _range=elements(TXh->mesh()),
                                _backend=backend
                                );
    opI->matPtr()->printMatlab("opIM.m");


    auto u = vf::project( _space=Xh,
                          _range=elements(mesh),
                          _expr=cos( pi*Px() )*sin( pi*Py() ) );

    auto tu = TXh->element();
    opI->apply( u,tu );

    auto error = integrate( _range=elements( TXh->mesh() ),
                            _expr=(idv( u )-idv( tu) )*(idv( u )-idv( tu) ) ).evaluate()( 0,0 );

    std::cout<<"error= "<< error <<"\n";

} // Test::run

int
main( int argc, char** argv )
{

    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=about(_name="test_hyper_interp",
                                  _author="Abdoulaye Samake",
                                  _email="abdoulaye.samake@imag.fr") );
    Application app;

    app.add( new Test<2,1>() );
    app.run();

}
