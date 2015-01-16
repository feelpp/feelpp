/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Abdoulaye Samake <abdoulaye.samake1@e.ujf-grenoble.fr>
       Date: 2013-11-15

  Copyright (C) 2013 Feel++ Consortium

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
   \author Abdoulaye Samake <abdoulaye.samake.@imag.fr>
   \date 2013-11-15
 */

#define BOOST_TEST_MODULE test_operatorinterpolation_hypercube
#include <testsuite/testsuite.hpp>

#include <feel/options.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/operatortrace.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>



namespace test_interpolation
{
using namespace Feel;

inline
AboutData
makeAbout()
{
    AboutData about( "Test_OperatorInterpolationHypercube" ,
                     "Test_OperatorInterpolationHypercube" ,
                     "0.1",
                     "test operator interpolation on hypercubes",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2013 Feel++ Consortium" );

    about.addAuthor( "Abdoulaye Samake", "developer", "asamake@gmail.com", "" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;
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
        backend( backend_type::build( soption( _name="backend" ) ) )
    {}

    void run();

private:

    backend_ptrtype backend;
}; // Test

template<int Dim, int Order>
void
Test<Dim,Order>::run()
{

    std::cout << "------------------------------------------------------------\n";
    std::cout << "Execute Test<" << Dim << "," << Order << ">\n";

    Environment::changeRepository( boost::format( "testsuite/feelinterpolation/%1%/%2%-%3%/P%4%/h_%5%/" )
                                   % this->about().appName()
                                   % convex_type::name()
                                   % Dim
                                   % Order
                                   % doption(_name="gmsh.hsize") );

    auto mesh = createGMSHMesh( _mesh=new mesh_type,
                                _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                                _desc=domain( _name=( boost::format( "hypercube-%1%" ) % Dim ).str(),
                                              _addmidpoint=false,
                                              _usenames=false,
                                              _shape="hypercube",
                                              _dim=Dim,
                                              _h=doption(_name="gmsh.hsize") ,
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
                          _range=markedfaces(mesh,"NORTH"),
                          _expr=cst(1.) );
                          //_expr=cos( pi*Px() )*sin( pi*Py() ) );

    auto tu = TXh->element();
    opI->apply( u,tu );
    u.printMatlab( "u1.m" );
    tu.printMatlab( "tu1.m" );

    auto l2error = normL2( _range=elements( TXh->mesh() ), _expr=idv( u )-idv( tu) );

    BOOST_TEST_MESSAGE( "l2 error 1 = " << l2error );
    BOOST_CHECK_SMALL( l2error, 1e-13 );

    u = vf::project( _space=Xh,
                     _range=markedfaces(mesh,"NORTH"),
                     _expr=sin( pi*Px() )*cos( pi*Py() ) );
    opI->apply( u,tu );
    u.printMatlab( "u2.m" );
    tu.printMatlab( "tu2.m" );

    l2error = normL2( _range=elements( TXh->mesh() ), _expr=idv( u )-idv( tu) );

    BOOST_TEST_MESSAGE( "l2 error 2 = " << l2error );
    BOOST_CHECK_SMALL( l2error, 1e-13 );

} // Test::run

}
/**
 * main code
 */
FEELPP_ENVIRONMENT_WITH_OPTIONS( test_interpolation::makeAbout(),
                                 Feel::feel_options() )

BOOST_AUTO_TEST_SUITE( interp_operatorinterpolation )

BOOST_AUTO_TEST_CASE( interp_operatorinterpolation_2d_1d_geo1 )
{
    BOOST_TEST_MESSAGE( "interp_interpolation_hypercube_2d_1d_q1" );
    using namespace test_interpolation;

    Test<2,1> t;
    t.run();
    BOOST_TEST_MESSAGE( "interp_interpolation_hypercube_2d_1d_q1 done" );
}

BOOST_AUTO_TEST_SUITE_END()
