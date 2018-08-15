/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2012-04-11

  Copyright (C) 2012 Université Joseph Fourier (Grenoble I)

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
   \file test_stencil.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2012-04-11
 */

// make sure that the init_unit_test function is defined by UTF
//#define BOOST_TEST_MAIN
// give a name to the testsuite
#define BOOST_TEST_MODULE stencil
#include <testsuite.hpp>

#include <feel/feelcore/application.hpp>
#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/pchv.hpp>
#include <feel/feeldiscr/thch.hpp>


#include <feel/feelvf/vf.hpp>

#include <feel/feeldiscr/projector.hpp>



namespace Feel
{
/**
 * This routine returns the list of options using the
 * boost::program_options library. The data returned is typically used
 * as an argument of a Feel::Application subclass.
 *
 * \return the list of options
 */
inline
po::options_description
makeOptions()
{
    po::options_description teststenciloptions( "test stencil options" );
    teststenciloptions.add_options()
    ( "hsize", po::value<double>()->default_value( 0.1 ), "mesh size" )
    ;
    return teststenciloptions.add( Feel::feel_options() );
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_stencil" ,
                     "test_stencil" ,
                     "0.1",
                     "Test for stencil manager",
                     AboutData::License_GPL,
                     "Copyright (c) 2012 Universite Joseph Fourier"
                     "Copyright (c) 2012-2016 Feel++ Consortium"
                     );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}

class TestStencil
    :
public Application
{
    typedef Application super;

public:

    //! numerical type is double
    typedef double value_type;

    //! linear algebra backend factory
    typedef Backend<value_type> backend_type;
    //! linear algebra backend factory shared_ptr<> type
    typedef typename std::shared_ptr<backend_type> backend_ptrtype ;

    //! geometry entities type composing the mesh, here Simplex in Dimension Dim of Order G_order
    typedef Simplex<2,1> convex_type;
    //! mesh type
    typedef Mesh<convex_type> mesh_type;
    //! mesh shared_ptr<> type
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    //! the basis type of our approximation space
    typedef bases<Lagrange<1> > basis_type;
    //! the approximation function space type
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef std::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;

    /**
     * Constructor
     */
    TestStencil()
        :
        super(),
        M_backend( backend_type::build( soption( _name="backend" ) ) ),
        meshSize( this->vm()["hsize"].as<double>() )
    {
        std::cout << "[TestStencil]\n";

        this->changeRepository( boost::format( "%1%/h_%2%/" )
                                % this->about().appName()
                                % this->vm()["hsize"].as<double>()
                              );



    }
    void
    operator()()
        {
            auto mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _desc = domain(_name="square",_h=meshSize,_shape="hypercube" ) );
            auto Xh = space_type::New( mesh );
            stencilManagerPrint();
            {
                auto M = M_backend->newMatrix( Xh, Xh );
                stencilManagerPrint();
                {
                    auto P = M_backend->newMatrix( Xh, Xh );
                    stencilManagerPrint();
                }
            }

            BOOST_CHECK_EQUAL( StencilManager::instance().size(), 1 );
            std::cout << "run garbage collector\n";
            stencilManagerGarbageCollect();
            BOOST_CHECK_EQUAL( StencilManager::instance().size(), 0 );

            auto M = M_backend->newMatrix( Xh, Xh );
            BOOST_CHECK_EQUAL( StencilManager::instance().size(), 1 );

            auto P = M_backend->newMatrix( Xh, Xh );
            BOOST_CHECK_EQUAL( StencilManager::instance().size(), 1 );

            stencilManagerPrint();

        }
private:
    //! linear algebra backend
    backend_ptrtype M_backend;

    //! mesh characteristic size
    double meshSize;

}; //TestStencil

struct A
{
    static const int nSpaces = 1;

    single_space_t f( single_space_t )
        {
            return single_space_t{};
        }
};
struct B
{
    static const int nSpaces = 3;

    multiple_space_t f( multiple_space_t )
        {
            return multiple_space_t{};
        }
};
}
#if USE_BOOST_TEST
FEELPP_ENVIRONMENT_WITH_OPTIONS( Feel::makeAbout(), Feel::makeOptions() )

BOOST_AUTO_TEST_SUITE( space )
BOOST_AUTO_TEST_CASE( test_stencil1 )
{
    BOOST_TEST_MESSAGE( "test_stencil1\n" );
    using namespace Feel;
    using namespace hana::literals;
    BOOST_HANA_CONSTANT_CHECK( single_space_t{} == 1_c);
    BOOST_HANA_CONSTANT_CHECK( multiple_space_t{} == 2_c);
    BOOST_HANA_CONSTANT_CHECK( single_spaces_t{} == 1_c);
    BOOST_HANA_CONSTANT_CHECK( multiple_spaces_t{} ==  2_c);
}
BOOST_AUTO_TEST_CASE( test_stencil2 )
{
    BOOST_TEST_MESSAGE( "test_stencil2\n" );
    using namespace Feel;
    using namespace hana::literals;

    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>> );
    auto Xh = Pch<1>( mesh );
    auto Yh = THch<1>( mesh );

    BOOST_HANA_CONSTANT_CHECK( type_space_t( Xh ) == single_space_c );
    BOOST_HANA_CONSTANT_CHECK( type_space_t( Yh ) == multiple_space_c );
    A a; B b;
    BOOST_HANA_CONSTANT_CHECK( type_space_t( a ) == single_space_c );
    BOOST_HANA_CONSTANT_CHECK( type_space_t( a ) == single_space_t() );
    BOOST_HANA_CONSTANT_CHECK( type_space_t( a ) == 1_c );
    BOOST_HANA_CONSTANT_CHECK( type_space_t( b ) == 2_c );
    std::cout << "type_space_t(a)=" << type_space_t(a) << "," << single_space_t{} << std::endl;
    BOOST_HANA_CONSTANT_CHECK( a.f( single_space_c ) == single_space_t{} );
    BOOST_HANA_CONSTANT_CHECK( single_space_c == single_space_t{} );
    BOOST_HANA_CONSTANT_CHECK(hana::integral_c<int, 2> == hana::int_c<2>);
    //static_assert(std::is_same<decltype(single_space_c), decltype(type_space_t(a))>::value,
    // "decltype(single_space_c) != single_space_t");
#if 0

    static_assert(std::is_same<decltype(type_space_t(a)), single_space_t>::value,
                  "decltype(type_space_t(a)) != single_space_t");
    #endif
    BOOST_HANA_CONSTANT_CHECK( a.f( type_space_t( a ) ) == single_space_t{} );
    BOOST_HANA_CONSTANT_CHECK( b.f( type_space_t( b ) ) == multiple_space_t{} );
    BOOST_TEST_MESSAGE( "test_stencil2 done\n" );
}

BOOST_AUTO_TEST_CASE( test_stencil )
{
    BOOST_TEST_MESSAGE( "test_stencil" );
    Feel::TestStencil t;

    t();

    Feel::stencilManagerPrint();

    std::cout << "run garbage collector outside Feel\n";
    Feel::stencilManagerGarbageCollect();
    BOOST_CHECK_EQUAL( Feel::StencilManager::instance().size(), 0 );

    Feel::stencilManagerGarbageCollect();
    BOOST_CHECK_EQUAL( Feel::StencilManager::instance().size(), 0 );

    BOOST_TEST_MESSAGE( "test_stencil" );
}

BOOST_AUTO_TEST_CASE( test_stencil_zero )
{
    using namespace Feel;
    using mesh_t = Mesh<Simplex<2>>;
    auto mesh = loadMesh( new mesh_t );

    auto Vh = Pchv<1>( mesh, true );
    auto Wh = Pch<1>( mesh, true );
    BlocksBaseGraphCSR zero_graph(2,2);
    zero_graph(0,0) = stencil( _test=Vh,_trial=Vh, _diag_is_nonzero=false, _close=false,_pattern=(size_type)Pattern::ZERO)->graph();
    zero_graph(1,0) = stencil( _test=Wh,_trial=Vh, _diag_is_nonzero=false, _close=false,_pattern=(size_type)Pattern::ZERO)->graph();
    zero_graph(0,1) = stencil( _test=Vh,_trial=Wh, _diag_is_nonzero=false, _close=false,_pattern=(size_type)Pattern::ZERO)->graph();
    zero_graph(1,1) = stencil( _test=Wh,_trial=Wh, _diag_is_nonzero=false, _close=false,_pattern=(size_type)Pattern::ZERO)->graph();

    auto A = backend()->newBlockMatrix(_block=zero_graph);
}

BOOST_AUTO_TEST_SUITE_END()
#else

int
main( int argc, char* argv[] )
{

    Feel::TestStencil app_stencil( argc, argv, Feel::makeAbout(), Feel::makeOptions() );

    // app_stencil.tangent_operators();
    app_stencil.shape_functions();
    // app_stencil.matrix_assembly();
    // app_stencil.example_problem();
}

#endif
