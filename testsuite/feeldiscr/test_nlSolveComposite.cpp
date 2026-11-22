/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Cecile Daversin <daversin@math.unistra.fr>
       Date: 2014-05-06

  Copyright (C) 2011 UJF
  Copyright (C) 2011 CNRS

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file test_nlSolveComposite.cpp
   \author Cecile Daversin <daversin@math.unistra.fr>
   \date 2014-01-29
 */
#define USE_BOOST_TEST 1

// make sure that the init_unit_test function is defined by UTF
//#define BOOST_TEST_MAIN
// give a name to the testsuite
#define BOOST_TEST_MODULE nlSolve PkPk
// disable the main function creation, use our own
//#define BOOST_TEST_NO_MAIN

#include <feel/feelcore/testsuite.hpp>
/** include predefined feel command line options */
#include <feel/options.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feeldiscr/mesh.hpp>
// #include <feel/feelfilters/creategmshmesh.hpp>
// #include <feel/feelfilters/loadgmshmesh.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/geo.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/product.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/blockforms.hpp>
#include <feel/feelvf/vf.hpp>
#include <cmath>
#include <algorithm>
#include <array>
#include <tuple>
#include <boost/test/data/test_case.hpp>

using namespace Feel;
namespace bdata = boost::unit_test::data;

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
    po::options_description testnlSolveoptions( "test_nlSolve options" );
    testnlSolveoptions.add_options()
    ( "hsize", po::value<double>()->default_value( 0.1 ), "mesh size" )
    ( "xmin", po::value<double>()->default_value( -1 ), "xmin of the reference element" )
    ( "ymin", po::value<double>()->default_value( -1 ), "ymin of the reference element" )
    ( "zmin", po::value<double>()->default_value( -1 ), "zmin of the reference element" )
    ( "sigma0", po::value<double>()->default_value( 1e+4 ), "mesh size" )
    ( "k0", po::value<double>()->default_value( 30 ), "mesh size" )
    ( "alpha", po::value<double>()->default_value( 4e-3 ), "mesh size" )
    ( "h_ech", po::value<double>()->default_value( 1000 ), "mesh size" )
    ( "Tw", po::value<double>()->default_value( 300 ), "mesh size" )
    ;
    return testnlSolveoptions.add( Feel::feel_options() );
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_nlSolveComposite" ,
                     "test_nlSolveComposite" ,
                     "0.1",
                     "Test for nlSolve with composite space",
                     AboutData::License_GPL,
                     "Copyright (c) 2009 Universite Joseph Fourier" );
    about.addAuthor( "Cecile Daversin", "developer", "daversin@math.unistra.fr", "" );
    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}

template<int Dim, int OrderV, int OrderT>
class TestNLSolveComposite
{
public:
    static constexpr int IntOrder_k = 2;
    static constexpr int IntOrder_dk = 1;
    static constexpr int maxNewtonIter = 20;
    static constexpr double newtonTol = 1e-12;

    using value_type = double;
    using backend_type = Backend<value_type>;
    using backend_ptrtype = std::shared_ptr<backend_type>;

    using convex_type = Simplex<Dim,1>;
    using mesh_type = Mesh<convex_type>;
    using mesh_ptrtype = std::shared_ptr<mesh_type>;

    using V_space_type = Pch_type<mesh_type, OrderV, value_type, PointSetFekete>;
    using T_space_type = Pch_type<mesh_type, OrderT, value_type, PointSetFekete>;
    using V_space_ptrtype = Pch_ptrtype<mesh_type, OrderV, value_type, PointSetFekete>;
    using T_space_ptrtype = Pch_ptrtype<mesh_type, OrderT, value_type, PointSetFekete>;

    TestNLSolveComposite()
        :
        M_backend( backend(_rebuild=true) ),
        meshSize( doption(_name="hsize") ),
        sigma0( doption(_name="sigma0") ),
        k0( doption(_name="k0") ),
        alpha( doption(_name="alpha") ),
        h( doption(_name="h_ech") ),
        Tw( doption(_name="Tw") )
    {
        Environment::changeRepository( _directory = boost::format( "testsuite/feeldiscr/test_nlSolveComposite/P%1%P%2%/h_%3%/" )
                                                    % OrderV % OrderT % meshSize );

        mesh = loadMesh( _mesh=new mesh_type );
        M_Vh = Pch<OrderV, value_type, PointSetFekete, mesh_type>( mesh );
        M_Th = Pch<OrderT, value_type, PointSetFekete, mesh_type>( mesh );
    }

    void run();

private:
    static uint16_type qOrder( int order )
    {
        return static_cast<uint16_type>( std::max( 1, order ) );
    }

    backend_ptrtype M_backend;
    mesh_ptrtype mesh;
    V_space_ptrtype M_Vh;
    T_space_ptrtype M_Th;
    double meshSize;
    double sigma0, k0, alpha, h, Tw;
};

template<int Dim, int OrderV, int OrderT>
void
TestNLSolveComposite<Dim, OrderV, OrderT>::run()
{
    auto ps = product( M_Vh, M_Th );
    auto U = ps.element();

    U(0_c).on( _range=elements(mesh), _expr=cst(1.) );
    U(1_c).on( _range=elements(mesh), _expr=cst(293.) );

    bool converged = false;
    for ( int kiter = 0; kiter < maxNewtonIter; ++kiter )
    {
        auto V = U(0_c);
        auto T = U(1_c);
        auto v = ps[0_c]->element();
        auto t = ps[1_c]->element();

        auto J = blockform2( ps, solve::strategy::monolithic, M_backend );
        auto R = blockform1( ps, solve::strategy::monolithic, M_backend );
        auto R0 = R(0_c);
        auto R1 = R(1_c);
        auto J00 = J(0_c,0_c);
        auto J01 = J(0_c,1_c);
        auto J10 = J(1_c,0_c);
        auto J11 = J(1_c,1_c);

        auto constexpr penaldir = 50;
        auto constexpr pV = OrderV;
        auto constexpr pT = OrderT;

        auto T0 = cst(293.);
        auto sigma = cst(sigma0)/( cst(1.) + cst(alpha)*(idv(T)-T0) );
        auto kval = (cst(k0)/(sigma0*T0))*sigma*idv(T);
        auto sigma_prime = -cst(alpha)*cst(sigma0)/(( cst(1.)+cst(alpha)*(idv(T)-T0))*(cst(1.)+cst(alpha)*(idv(T)-T0)) );
        auto k_prime = (cst(k0)/(cst(sigma0)*T0))*(sigma_prime*idv(T)+sigma);

        auto Dirichlet_cst = 1+Py();
        auto Robin_cst = -h*Tw;
        auto Robin_coeff = h;

        // RHS = -R(U)
        R0 += integrate( _range=elements(mesh),
                         _expr=-val(sigma)*grad(v)*trans(gradv(V)),
                         _quad=_Q( qOrder( IntOrder_k + (pV-1)+(pV-1) ) ) );
        R1 += integrate( _range=elements(mesh),
                         _expr=-val(kval)*grad(t)*trans(gradv(T)),
                         _quad=_Q( qOrder( IntOrder_k + (pT-1)+(pT-1) ) ) );
        R1 += integrate( _range=elements(mesh),
                         _expr=-val(-sigma*gradv(V)*trans(gradv(V)))*id(t),
                         _quad=_Q( qOrder( IntOrder_k + (pV-1)+(pV-1)+pT ) ) );

        R0 += integrate( _range=boundaryfaces(mesh),
                         _expr=-val(-sigma*gradv(V)*vf::N())*id(v),
                         _quad=_Q( qOrder( IntOrder_k + (pV-1)+pV ) ) );
        R0 += integrate( _range=boundaryfaces(mesh),
                         _expr=-val(sigma*penaldir*idv(V)/hFace())*id(v),
                         _quad=_Q( qOrder( IntOrder_k + pV+pV ) ) );
        R0 += integrate( _range=boundaryfaces(mesh),
                         _expr=-grad(v)*val(-sigma*vf::N()*idv(V)),
                         _quad=_Q( qOrder( IntOrder_k + (pV-1)+pV ) ) );
        R0 += integrate( _range=boundaryfaces(mesh),
                         _expr=-val(-sigma*penaldir/hFace())*Dirichlet_cst*id(v),
                         _quad=_Q( qOrder( IntOrder_k + pV ) ) );
        R0 += integrate( _range=boundaryfaces(mesh),
                         _expr=-grad(v)*val(sigma*vf::N())*Dirichlet_cst,
                         _quad=_Q( qOrder( IntOrder_k + (pV-1) ) ) );
        R1 += integrate( _range=boundaryfaces(mesh),
                         _expr=-(Robin_coeff*idv(T)*id(t) + cst(Robin_cst)*id(t)) );

        // Jacobian blocks
        J00 += integrate( _range=elements(mesh),
                          _expr=val(sigma)*grad(v)*trans(gradt(v)),
                          _quad=_Q( qOrder( IntOrder_k + (pV-1)+(pV-1) ) ) );
        J01 += integrate( _range=elements(mesh),
                          _expr=val(sigma_prime*gradv(V))*trans(grad(v))*idt(t),
                          _quad=_Q( qOrder( IntOrder_dk + (pV-1)+(pV-1)+pT ) ) );
        J10 += integrate( _range=elements(mesh),
                          _expr=val(-2*sigma*gradv(V))*trans(gradt(v))*id(t),
                          _quad=_Q( qOrder( IntOrder_k + pV + (pV-1)*(pT-1) ) ) );
        J11 += integrate( _range=elements(mesh),
                          _expr=val(k_prime*gradv(T))*trans(grad(t))*idt(t),
                          _quad=_Q( qOrder( IntOrder_dk + (pT-1)+(pT-1)+pT ) ) );
        J11 += integrate( _range=elements(mesh),
                          _expr=val(kval)*grad(t)*trans(gradt(t)),
                          _quad=_Q( qOrder( IntOrder_k + (pT-1)+(pT-1) ) ) );
        J11 += integrate( _range=elements(mesh),
                          _expr=val(-sigma_prime*gradv(V)*trans(gradv(V)))*id(t)*idt(t),
                          _quad=_Q( qOrder( IntOrder_dk + (pV-1)+(pV-1)+pT+pT ) ) );

        J00 += integrate( _range=boundaryfaces(mesh),
                          _expr=gradt(v)*val(-sigma*vf::N())*id(v),
                          _quad=_Q( qOrder( IntOrder_k + (pV-1)+pV ) ) );
        J00 += integrate( _range=boundaryfaces(mesh),
                          _expr=sigma*penaldir*id(v)*idt(v)/hFace() );
        J00 += integrate( _range=boundaryfaces(mesh),
                          _expr=grad(v)*val(-sigma*vf::N())*idt(v),
                          _quad=_Q( qOrder( IntOrder_k + (pV-1)+pV ) ) );
        J01 += integrate( _range=boundaryfaces(mesh),
                          _expr=val(-sigma_prime*gradv(V)*vf::N())*id(v)*idt(t),
                          _quad=_Q( qOrder( IntOrder_dk + (pV-1)+pV+pT ) ) );
        J01 += integrate( _range=boundaryfaces(mesh),
                          _expr=grad(v)*val(-sigma_prime*vf::N()*idv(V))*idt(t),
                          _quad=_Q( qOrder( IntOrder_dk + (pV-1)+pV+pT ) ) );
        J01 += integrate( _range=boundaryfaces(mesh),
                          _expr=grad(v)*idt(t)*val(sigma_prime*vf::N()*Dirichlet_cst),
                          _quad=_Q( qOrder( IntOrder_dk + (pV-1)+pT ) ) );
        J01 += integrate( _range=boundaryfaces(mesh),
                          _expr=sigma_prime*penaldir*idv(V)*id(v)*idt(t)/hFace() );
        J01 += integrate( _range=boundaryfaces(mesh),
                          _expr=-sigma_prime*penaldir*Dirichlet_cst*id(v)*idt(t)/hFace() );
        J11 += integrate( _range=boundaryfaces(mesh),
                          _expr=Robin_coeff*id(t)*idt(t) );

        J.close();
        R.close();

        auto dU = ps.element();
        J.solve( _solution=dU, _rhs=R );
        U(0_c) += dU(0_c);
        U(1_c) += dU(1_c);

        auto const dVnorm = integrate( _range=elements(mesh), _expr=idv(dU(0_c))*idv(dU(0_c)) ).evaluate()(0,0);
        auto const dTnorm = integrate( _range=elements(mesh), _expr=idv(dU(1_c))*idv(dU(1_c)) ).evaluate()(0,0);
        auto const dUnorm = std::sqrt( dVnorm + dTnorm );
        if ( dUnorm < newtonTol )
        {
            converged = true;
            break;
        }
    }

    BOOST_CHECK_MESSAGE( converged, "Newton iterations did not converge" );

    auto const T_mean_int = integrate( _range=elements(mesh), _expr=idv(U(1_c)) ).evaluate()(0,0);
    auto const area = integrate( _range=elements(mesh), _expr=cst(1.) ).evaluate()(0,0);
    auto const T_mean = T_mean_int/area;

    if ( Environment::worldComm().globalRank() == Environment::worldComm().masterRank() )
        std::cout << "[P" << OrderV << "-P" << OrderT << "] Tmean = " << T_mean << std::endl;

    BOOST_CHECK_CLOSE( T_mean, 344, 2e-1 );

    auto e = exporter( _mesh=mesh, _name=( boost::format( "test_nlSolveComposite_P%1%P%2%" ) % OrderV % OrderT ).str() );
    if ( e->doExport() )
    {
        e->step(0)->setMesh( mesh );
        e->step(0)->add( "Potential", U(0_c) );
        e->step(0)->add( "Temperature", U(1_c) );
        e->save();
    }
}

#if USE_BOOST_TEST

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )
BOOST_AUTO_TEST_SUITE( NLSOLVE_COMPOSITE )

template<int OrderV, int OrderT>
inline void runNLSolveCase2D()
{
    BOOST_TEST_CONTEXT( "electro-thermal non linear model - 2D - P" << OrderV << "P" << OrderT )
    {
        TestNLSolveComposite<2, OrderV, OrderT> app_testNLSolve;
        app_testNLSolve.run();
    }
}

inline void runNLSolveCase2D( int orderV, int orderT )
{
    if ( orderV == 1 && orderT == 2 )
        runNLSolveCase2D<1,2>();
    else if ( orderV == 2 && orderT == 1 )
        runNLSolveCase2D<2,1>();
    else if ( orderV == 1 && orderT == 1 )
        runNLSolveCase2D<1,1>();
    else if ( orderV == 2 && orderT == 2 )
        runNLSolveCase2D<2,2>();
    else
        BOOST_FAIL( "unsupported (orderV, orderT) combination in nlSolve_2D dataset" );
}

BOOST_DATA_TEST_CASE( nlSolve_2D,
                      bdata::make( std::array<std::tuple<int,int>,4>{ std::tuple{1,2},
                                                                       std::tuple{2,1},
                                                                       std::tuple{1,1},
                                                                       std::tuple{2,2} } ),
                      orderV, orderT )
{
    runNLSolveCase2D( orderV, orderT );
}

BOOST_AUTO_TEST_SUITE_END()
#endif
