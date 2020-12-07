/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2012-04-27

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
   \file test_backend.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2012-04-27
 */
#define USE_BOOST_TEST 1
#define BOOST_TEST_MODULE backend testsuite

#include <feel/feelcore/testsuite.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/application.hpp>
#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
//#include <feel/feeldiscr/stencil.hpp>
#include <feel/feelvf/vf.hpp>

#define FEELAPP()                                                       \
    Feel::Application app;                                            \
    if ( app.vm().count( "help" ) )                                     \
    {                                                                   \
        std::cout << app.optionsDescription() << "\n";                  \
    }



namespace Feel
{
inline
po::options_description
makeOptions()
{
    po::options_description simgetoptions( "test_backend options" );
    simgetoptions.add_options()
    ( "hsize", po::value<double>()->default_value( 0.5 ), "mesh size" )
    ;
    return simgetoptions.add( Feel::feel_options() )
                        .add( backend_options("test1"))
                        .add( backend_options("ctest1"));
}


inline
AboutData
makeAbout()
{
    AboutData about( "test_backend" ,
                     "test_backend" ,
                     "0.1",
                     "Backend tests",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2012 Feel++ Consortium" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}

/**
 *
 */
class sim : public Simget
{
public:
    sim( )
        :
        Simget(),
        meshSize( this->vm()["hsize"].as<double>() )
        {
        }
    void run() override

        {
            BOOST_CHECK( mpi::environment::initialized() );
            BOOST_TEST_MESSAGE( "mpi ok" );
            BOOST_CHECK_EQUAL( Feel::detail::BackendManager<double>::instance().size(), 0 );
            BOOST_CHECK_EQUAL( Feel::detail::BackendManager<std::complex<double>>::instance().size(), 0 );
            BOOST_TEST_MESSAGE( "backend manager empty" );
            auto b1 = backend(_name="test1");
            auto c1 = cbackend(_name="ctest1");
#if 0
            BOOST_TEST_MESSAGE( "creating backend" );
            BOOST_CHECK_EQUAL( Feel::detail::BackendManager::instance().size(), 1 );
            BOOST_CHECK_EQUAL( Feel::detail::BackendManager::instance().begin()->first.first, BACKEND_PETSC );
            BOOST_CHECK_EQUAL( Feel::detail::BackendManager::instance().begin()->first.second, "test1" );
            BOOST_FOREACH( auto b, Feel::detail::BackendManager::instance() )
            {
                BOOST_TEST_MESSAGE( "[test_backend] backend name =" << b.first.second << " kind =" << b.first.first );;
            }
            auto v = b->newVector( 10, 10 );
#endif
        }
private:
    double meshSize;
};

namespace detail
{
template <typename T>
double myLocalProcessSum( Feel::Vector<T> const& vec )
{
    double res = 0;
    for ( size_type k=0;k<vec.map().nLocalDofWithGhost();++k )
        res += vec( k );
    return res;
}
}

} // Feel

using namespace Feel;

FEELPP_ENVIRONMENT_WITH_OPTIONS( Feel::makeAbout(), Feel::makeOptions() )

BOOST_AUTO_TEST_SUITE( backendsuite )
BOOST_AUTO_TEST_CASE( test_backend1 )
{

    BOOST_TEST_MESSAGE( "test_backend1" );
    BOOST_CHECK( Environment::initialized() );
    //BOOST_CHECK( mpi::environment::initialized() );
    BOOST_TEST_MESSAGE( "initializing the Application" );

    FEELAPP();
    BOOST_CHECK( mpi::environment::initialized() );
    BOOST_TEST_MESSAGE( "adding simget" );
    app.add( new sim );
    BOOST_TEST_MESSAGE( "run simget" );
    app.run();

    BOOST_TEST_MESSAGE( "test_backend1 done" );
}
BOOST_AUTO_TEST_CASE( test_backend_petsc )
{
    double tolCheck = 1e-9;
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>);
    double meshMeasure = mesh->measure();
    auto Vh1 = Pch<2>( mesh );
    auto u1 = Vh1->element();
    size_type nDofVh1 = Vh1->nDof();
    size_type nLocalDofWithGhostVh1 = Vh1->nLocalDofWithGhost();
    size_type nLocalDofWithoutGhostVh1 = Vh1->nLocalDofWithoutGhost();

    auto backendPetsc = backend(_kind="petsc");
    auto vec_petsc1x = backendPetsc->newVector( Vh1 );
    auto vec_petsc1y = backendPetsc->newVector( Vh1 );
    auto vec_petsc1z = backendPetsc->newVector( Vh1 );
    auto mat_petsc1 = backendPetsc->newMatrix( _test=Vh1, _trial=Vh1 );
    auto mat_petsc1y = backendPetsc->newMatrix( _test=Vh1, _trial=Vh1 );

    // assembly matrix, vector
    form2(_test=Vh1,_trial=Vh1,_matrix=mat_petsc1 ) =
        integrate(_range=elements(mesh),_expr=idt(u1)*id(u1) );
    mat_petsc1->close();
    form1(_test=Vh1,_vector=vec_petsc1z ) =
        integrate(_range=elements(mesh),_expr=id(u1) );
    vec_petsc1z->close();

    vec_petsc1x->setConstant( 3. );
    u1.setConstant( 2. );

    // pre-check
    BOOST_CHECK( vec_petsc1x->sum() == 3*nDofVh1 );
    //BOOST_CHECK_CLOSE( mat_petsc1->energy(u1,u1) , 3*meshMeasure , tolCheck );
    BOOST_CHECK_CLOSE( mat_petsc1->energy(vec_petsc1x,vec_petsc1x) , 3*3*meshMeasure , tolCheck );
    // prod
    backendPetsc->prod( mat_petsc1, vec_petsc1x, vec_petsc1y );
    BOOST_CHECK_CLOSE( vec_petsc1y->sum() , 3*vec_petsc1z->sum(), tolCheck );
    // diag
    mat_petsc1y->zero();
    for ( size_type k=0;k<nLocalDofWithoutGhostVh1;++k )
        mat_petsc1y->set( k,k, 6. );
    mat_petsc1y->close();
    backendPetsc->diag( *mat_petsc1y, *vec_petsc1x );
    BOOST_CHECK_CLOSE( vec_petsc1x->sum() , 6*nDofVh1, tolCheck );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*vec_petsc1x), 6*nLocalDofWithGhostVh1, tolCheck );
    // diag with clear
    vec_petsc1x->clear();
    backendPetsc->diag( *mat_petsc1y, *vec_petsc1x );
    BOOST_CHECK_CLOSE( vec_petsc1x->sum() , 6*nDofVh1, tolCheck );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(*vec_petsc1x), 6*nLocalDofWithGhostVh1, tolCheck );
    // diag with ublas vector
    backendPetsc->diag( *mat_petsc1y, u1 );
    BOOST_CHECK_CLOSE( u1.sum() , 6*nDofVh1, tolCheck );
    BOOST_CHECK_CLOSE( Feel::detail::myLocalProcessSum(u1), 6*nLocalDofWithGhostVh1, tolCheck );
    // set diag
    vec_petsc1x->setConstant( 4. );
    vec_petsc1y->setConstant( 5. );
    backendPetsc->diag( *vec_petsc1x, *mat_petsc1y );
    BOOST_CHECK_CLOSE( mat_petsc1y->energy(vec_petsc1y,vec_petsc1y) , 5*4*5*nDofVh1 , tolCheck );


}
BOOST_AUTO_TEST_SUITE_END()

