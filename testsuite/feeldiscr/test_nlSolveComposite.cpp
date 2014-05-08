/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

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

#include <testsuite/testsuite.hpp>
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
#include <feel/feeldiscr/region.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;

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
    :
public Application
{
    typedef Application super;

public:
    static const int IntOrder_k = 2;
    static const int IntOrder_dk = 1;

    //! numerical type is double
    typedef double value_type;
    typedef TestNLSolveComposite<Dim, OrderV, OrderT> self_type;

    //! linear algebra backend factory
    typedef Backend<value_type> backend_type;
    //! linear algebra backend factory shared_ptr<> type
    typedef typename boost::shared_ptr<backend_type> backend_ptrtype ;

    //! sparse matrix type associated with backend
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    //! geometry entities type composing the mesh, here Simplex in Dimension Dim of Order G_order
    typedef Simplex<Dim,1> convex_type;
    //! mesh type
    typedef Mesh<convex_type> mesh_type;
    //! mesh shared_ptr<> type
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef Lagrange<OrderV,Scalar> V_single_basis_type;
    typedef Lagrange<OrderT,Scalar> T_single_basis_type;
    typedef bases<V_single_basis_type, T_single_basis_type> basis_type;

    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    typedef typename space_type::element_type element_type;
    typedef boost::shared_ptr<element_type> element_ptrtype;

    //! the exporter factory type
    typedef Exporter<mesh_type> export_type;
    //! the exporter factory (shared_ptr<> type)
    typedef boost::shared_ptr<export_type> export_ptrtype;

    /**
     * Constructor
     */
    TestNLSolveComposite()
        :
        super(),
        M_backend( backend(_rebuild=true) ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        sigma0(this->vm()["sigma0"].template as<double>()),
        k0(this->vm()["k0"].template as<double>()),
        alpha(this->vm()["alpha"].template as<double>()),
        h(this->vm()["h_ech"].template as<double>()),
        Tw(this->vm()["Tw"].template as<double>())
    {
        this->changeRepository( boost::format( "%1%/h_%2%/" )
                                % this->about().appName()
                                % this->vm()["hsize"].template as<double>()
                                );

        mesh = loadMesh(_mesh = new mesh_type);

        M_backend->nlSolver()->jacobian = boost::bind( &self_type::updateJacobian,
                                                       boost::ref( *this ), _1, _2 );

        M_backend->nlSolver()->residual = boost::bind( &self_type::updateResidual,
                                                       boost::ref( *this ), _1, _2 );
    }

    /**
     * run the application
     */
    void updateResidual(const vector_ptrtype& X, vector_ptrtype& R);
    void updateJacobian(const vector_ptrtype& X, sparse_matrix_ptrtype& J);
    void newtonSolve(element_type& sol);
    void run();

private:
    //! linear algebra backend
    backend_ptrtype M_backend;

    //! mesh characteristic size
    mesh_ptrtype mesh;
    double meshSize;

    //! Weak Dirichlet treatment
    bool weakdir;
    double penaldir;

    double sigma0,k0,alpha,h,Tw;

}; //TestNLSolvePkPk

template<int Dim, int OrderV, int OrderT>
void
TestNLSolveComposite<Dim, OrderV, OrderT>::updateResidual(const vector_ptrtype& X, vector_ptrtype& R)
{
    space_ptrtype M_Xh = space_type::New(mesh);

    auto U = M_Xh->element(); //trial
    auto V_test = M_Xh->element(); //test

    U=*X;

    auto V = U.template element<0>(); //potential
    auto phi_V = V_test.template element<0>();

    auto T = U.template element<1>(); //temperature
    auto phi_T = V_test.template element<1>();

    R->zero();
    auto penaldir = 50;

    auto T0 = cst(293.);
    auto sigma = cst(sigma0)/( cst(1.) + cst(alpha)*(idv(T) - T0)  );
    auto k = (cst(k0)/(sigma0*T0))*sigma*idv(T);

    // Comp 1 : F1(V,T) = variationnal formulation for V
    form1( _test=M_Xh, _vector=R ) += integrate( _range=elements(mesh),
                                                 _expr=val(sigma)*grad(phi_V)*trans(gradv(V)),
                                                 _quad=_Q<IntOrder_k+(OrderV-1)+(OrderV-1)>()
                                                 );

    // Comp 2 : F2(V,T) = variationnal formulation for T
    form1( _test=M_Xh, _vector=R ) += integrate( _range=elements(mesh),
                                                 _expr=val(k)*grad(phi_T)*trans(gradv(T)),
                                                 _quad=_Q<IntOrder_k+(OrderT-1)+(OrderT-1)>()
                                                 );
    form1( _test=M_Xh, _vector=R ) += integrate( _range=elements(mesh),
                                                 _expr=val(-sigma*gradv(V)*trans(gradv(V)))*id(phi_T),
                                                 _quad=_Q<IntOrder_k+(OrderV-1)+(OrderV-1)+OrderT>()
                                                 );

    //// ***  Boundary conditions *** ////
    // Comp 1 : F1(V,T) = variationnal formulation for V

    auto Dirichlet_cst = 1+Py();
    auto Robin_cst = -h*Tw;
    auto Robin_coeff = h;

    form1( _test=M_Xh, _vector=R ) += integrate( _range=boundaryfaces(mesh),
                                                 _expr=val(-sigma*gradv(V)*vf::N())*id(phi_V),
                                                 _quad=_Q<IntOrder_k+(OrderV-1)+OrderV>()
                                                 );
    form1( _test=M_Xh, _vector=R ) += integrate( _range=boundaryfaces(mesh),
                                                 _expr=val(sigma*penaldir*idv(V)/hFace())*id(phi_V),
                                                 _quad=_Q<IntOrder_k+OrderV+OrderV>()
                                                 );
    form1( _test=M_Xh, _vector=R ) += integrate( _range=boundaryfaces(mesh),
                                                 _expr=grad(phi_V)*val(-sigma*vf::N()*idv(V)),
                                                 _quad=_Q<IntOrder_k+(OrderV-1)+OrderV>()
                                                 );
    form1( _test=M_Xh, _vector=R ) += integrate( _range=boundaryfaces(mesh),
                                                 _expr=val(-sigma*penaldir/hFace())*Dirichlet_cst*id(phi_V),
                                                 _quad=_Q<IntOrder_k+OrderV>()
                                                 );
    form1( _test=M_Xh, _vector=R ) += integrate( _range=boundaryfaces(mesh),
                                                 _expr=grad(phi_V)*val(sigma*vf::N())*Dirichlet_cst,
                                                 _quad=_Q<IntOrder_k+(OrderV-1)>()
                                                 );

    // Comp 2 : F1(V,T) = variationnal formulation for T
    form1( _test=M_Xh, _vector=R ) += integrate( _range=boundaryfaces(mesh),
                                                 _expr=Robin_coeff*idv(T)*id(phi_T)
                                                 +cst(Robin_cst)*id(phi_T) );
}

template<int Dim, int OrderV, int OrderT>
void
TestNLSolveComposite<Dim, OrderV, OrderT>::updateJacobian(const vector_ptrtype& X, sparse_matrix_ptrtype& J)
{
    space_ptrtype M_Xh = space_type::New(mesh);

    auto U = M_Xh->element(); //trial
    auto V_test = M_Xh->element(); //test

    U=*X;

    auto V = U.template element<0>(); //potential
    auto phi_V = V_test.template element<0>();

    auto T = U.template element<1>(); //temperature
    auto phi_T = V_test.template element<1>();

    J->zero();
    auto penaldir = 50;

    auto T0 = cst(293.);
    auto sigma = cst(sigma0)/( cst(1.) + cst(alpha)*(idv(T) - T0)  );
    auto k = (cst(k0)/(sigma0*T0))*sigma*idv(T);
    auto sigma_prime = -cst(alpha)*cst(sigma0)/(( cst(1.)+cst(alpha)*(idv(T)-T0))*(cst(1.)+cst(alpha)*(idv(T)-T0)) );
    auto k_prime = (cst(k0)/(cst(sigma0)*T0))*(sigma_prime*idv(T) + sigma );

    //// Comp(0,0) : d(F1)/dV
    form2( _test=M_Xh, _trial=M_Xh, _matrix=J ) += integrate( _range=elements( mesh ),
                                                              _expr=val(sigma)*grad(phi_V)*trans(gradt(V)),
                                                              _quad=_Q<IntOrder_k+(OrderV-1)+(OrderV-1)>()
                                                              );
    //// Comp(0,1) : d(F1)/dT
    form2( _test=M_Xh, _trial=M_Xh, _matrix=J ) += integrate( _range=elements( mesh ),
                                                              _expr=val(sigma_prime*gradv(V))*trans(grad(phi_V))*idt(T),
                                                              _quad=_Q<IntOrder_dk+(OrderV-1)+(OrderV-1)+OrderT>()
                                                              );
    //// Comp(1,0) : d(F2)/dV
    form2( _test=M_Xh, _trial=M_Xh, _matrix=J ) += integrate( _range=elements( mesh ),
                                                              _expr=val(-2*sigma*gradv(V))*trans(gradt(V))*id(phi_T),
                                                              _quad=_Q<IntOrder_k+OrderV+(OrderV-1)*(OrderT-1)>()
                                                              );

    //// Comp(1,1) : d(F2)/dT
    form2( _test=M_Xh, _trial=M_Xh, _matrix=J ) += integrate( _range=elements( mesh ),
                                                              _expr=val(k_prime*gradv(T))*trans(grad(phi_T))*idt(T),
                                                              _quad=_Q<IntOrder_dk+(OrderT-1)+(OrderT-1)+OrderT>()
                                                              );
    form2( _test=M_Xh, _trial=M_Xh, _matrix=J ) += integrate( _range=elements( mesh ),
                                                              _expr=val(k)*grad(phi_T)*trans(gradt(T)),
                                                              _quad=_Q<IntOrder_k+(OrderT-1)+(OrderT-1)>()
                                                              );
    form2( _test=M_Xh, _trial=M_Xh, _matrix=J ) += integrate( _range=elements( mesh ),
                                                              _expr=val(- sigma_prime*gradv(V))
                                                              *trans(gradt(V))*id(phi_T)*idt(T),
                                                              _quad=_Q<IntOrder_dk+(OrderV-1)+(OrderV-1)+OrderT+OrderT>()
                                                              );
    auto Dirichlet_cst = 1+Py();
    auto Robin_cst = -h*Tw;
    auto Robin_coeff = h;

    //// Comp(0,0) : d(F1)/dV
    form2( _test=M_Xh, _trial=M_Xh, _matrix=J ) += integrate( _range=boundaryfaces(mesh),
                                                              _expr=gradt(V)*val(- sigma*vf::N())*id(phi_V),
                                                              _quad=_Q<IntOrder_k+(OrderV-1)+OrderV>()
                                                              );
    form2( _test=M_Xh, _trial=M_Xh, _matrix=J ) += integrate( _range=boundaryfaces(mesh),
                                                              _expr=sigma*penaldir*id(phi_V)*idt(V)/hFace() );
    form2( _test=M_Xh, _trial=M_Xh, _matrix=J ) += integrate( _range=boundaryfaces(mesh),
                                                              _expr=grad(phi_V)*val(- sigma*vf::N())*idt(V),
                                                              _quad=_Q<IntOrder_k+(OrderV-1)+OrderV>()
                                                              );
    //// Comp(0,1) : d(F1)/dT
    form2( _test=M_Xh, _trial=M_Xh, _matrix=J ) += integrate( _range=boundaryfaces(mesh),
                                                              _expr=val(- sigma_prime*gradv(V)*vf::N())
                                                              *id(phi_V)*idt(T),
                                                              _quad=_Q<IntOrder_dk+(OrderV-1)+OrderV+OrderT>()
                                                              );
    form2( _test=M_Xh, _trial=M_Xh, _matrix=J ) += integrate( _range=boundaryfaces(mesh),
                                                              _expr=grad(phi_V)*val(- sigma_prime*vf::N()*idv(V))
                                                              *idt(T),
                                                              _quad=_Q<IntOrder_dk+(OrderV-1)+OrderV+OrderT>()
                                                              );
    form2( _test=M_Xh, _trial=M_Xh, _matrix=J ) += integrate( _range=boundaryfaces(mesh),
                                                              _expr=grad(phi_V)*idt(T)
                                                              *val(sigma_prime*vf::N()*Dirichlet_cst),
                                                              _quad=_Q<IntOrder_dk+(OrderV-1)+OrderT>()
                                                              );
    form2( _test=M_Xh, _trial=M_Xh, _matrix=J ) += integrate( _range=boundaryfaces(mesh),
                                                              _expr = sigma_prime*penaldir*idv(V)*id(phi_V)*idt(T)/hFace() );
    form2( _test=M_Xh, _trial=M_Xh, _matrix=J ) += integrate( _range=boundaryfaces(mesh),
                                                              _expr = -sigma_prime*penaldir*Dirichlet_cst
                                                              *id(phi_V)*idt(T)/hFace() );

    //// Comp(1,1) : d(F2)/dT
    form2( _test=M_Xh, _trial=M_Xh, _matrix=J ) += integrate( _range=boundaryfaces(mesh),
                                                              _expr=Robin_coeff*id(phi_T)*idt(T) );
}

template<int Dim, int OrderV, int OrderT>
void
TestNLSolveComposite<Dim, OrderV, OrderT>::newtonSolve(element_type& sol)
{
    space_ptrtype M_Xh = sol.functionSpace();

    auto U_vec = M_backend->newVector(M_Xh);
    *U_vec = sol;

    auto J = M_backend->newMatrix(_test=M_Xh,_trial=M_Xh);
    auto R = M_backend->newVector(M_Xh);

    U_vec->close();

    M_backend->nlSolve(_jacobian=J, _solution=U_vec, _residual=R);

    sol = *U_vec;
}

template<int Dim, int OrderV, int OrderT>
void
TestNLSolveComposite<Dim, OrderV, OrderT>::run()
{
    auto Xh = space_type::New(mesh);
    auto VT = Xh->element();

    auto V = vf::project( VT.template element<0>().functionSpace(), elements(mesh), cst(1.) );
    auto T = vf::project( VT.template element<1>().functionSpace(), elements(mesh), cst(293.) );

    VT.template element<0>() = V;
    VT.template element<1>() = T;

    newtonSolve(VT);

    V = VT.template element<0>();
    T = VT.template element<1>();

    auto T_mean = integrate( elements(mesh), idv(T) ).evaluate()(0,0);
    auto area = integrate( elements(mesh), cst(1.) ).evaluate()(0,0);
    T_mean /= area;
    if( Environment::worldComm().globalRank() == Environment::worldComm().masterRank() )
        std::cout << "[P" << OrderV << "-P" << OrderT << "] Tmean = " << T_mean << std::endl;

    BOOST_CHECK_CLOSE( T_mean, 344, 2e-1 );

    TestNLSolveComposite::export_ptrtype exporter( TestNLSolveComposite::export_type::New( this->vm(),
                                                                                           (boost::format( "%1%" )
                                                                                            % this->about().appName() ).str() ) );

    if ( exporter->doExport() )
        {
            exporter->step(0)->setMesh( mesh );
            exporter->step(0)->add( "Potential", VT.template element<0>() );
            exporter->step(0)->add( "Temperature", VT.template element<1>() );
            exporter->save();
        }
}

#if USE_BOOST_TEST

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )
BOOST_AUTO_TEST_SUITE( NLSOLVE_COMPOSITE )

BOOST_AUTO_TEST_CASE( nlSolve_2D_P1P2)
{
    BOOST_TEST_MESSAGE( "*** electro-thermal non linear model - 2D - P1P2 ***" );

    TestNLSolveComposite<2,1,2> app_testNLSolve;
    app_testNLSolve.run();
}

BOOST_AUTO_TEST_CASE( nlSolve_2D_P2P1)
{
    BOOST_TEST_MESSAGE( "*** electro-thermal non linear model - 2D - P2P1 ***" );

    TestNLSolveComposite<2,2,1> app_testNLSolve;
    app_testNLSolve.run();
}

BOOST_AUTO_TEST_CASE( nlSolve_2D_P1P1)
{
    BOOST_TEST_MESSAGE( "*** electro-thermal non linear model - 2D - P1P1 ***" );

    TestNLSolveComposite<2,1,1> app_testNLSolve;
    app_testNLSolve.run();
}

BOOST_AUTO_TEST_CASE( nlSolve_2D_P2P2)
{
    BOOST_TEST_MESSAGE( "*** electro-thermal non linear model - 2D - P2P2 ***" );

    TestNLSolveComposite<2,2,2> app_testNLSolve;
    app_testNLSolve.run();
}

BOOST_AUTO_TEST_SUITE_END()
#endif
