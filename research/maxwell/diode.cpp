/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2011-06-01

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
   \file diode.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2011-06-01
 */
#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/projector.hpp>
#include <feel/feelvf/vf.hpp>

/** use Feel namespace */
using namespace Feel;
using namespace Feel::vf;

enum RKMethod
    {
        EULER_EXPLICIT = 0,
        EULER_MODIFIED,
        HEUN,
        RK4
    };

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
    po::options_description diodeoptions("Diode options");
    diodeoptions.add_options()
        ("hsize", po::value<double>()->default_value( 0.5 ), "mesh size")
        ("convex", Feel::po::value<std::string>()->default_value( "simplex" ), "shape of the convex used in the mesh (either simplex or hypercube)")
        ("penaldir", Feel::po::value<double>()->default_value( 20 ), "penalisation parameter for the weak boundary conditions")
        ("dt", Feel::po::value<double>()->default_value( 0.01 ), "timestep value")
        ("Tfinal", Feel::po::value<double>()->default_value( 1 ), "final time")
        ("rkmethod", Feel::po::value<int>()->default_value( 2 ), "rk method, 0=euler, 1=modified euler, 2=heun, 3=rk4")
        ;
    return diodeoptions.add( Feel::feel_options() );
}

/**
 * This routine defines some information about the application like
 * authors, version, or name of the application. The data returned is
 * typically used as an argument of a Feel::Application subclass.
 *
 * \return some data about the application.
 */
inline
AboutData
makeAbout()
{
    AboutData about( "diode" ,
                     "diode" ,
                     "0.1",
                     "2D Diode",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2011 Universite Joseph Fourier");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    about.addAuthor("Philippe Helluy", "developer", "helluy@math.unistra.fr", "");
    return about;

}

template<typename A, uint16_type i>
class BasisTag : public A
{
public:
    static const uint16_type TAG = i;

};

/**
 * diode geometry description
 */
gmsh_ptrtype diodegeo( double h, int Order, std::string const& convex );

/**
 * \class Diode
 *
 * Diode Solver using continuous approximation spaces
 * solve \f$ -\Delta u = f\f$ on \f$\Omega\f$ and \f$u= g\f$ on \f$\Gamma\f$
 *
 * \tparam Dim the geometric dimension of the problem (e.g. Dim=1, 2 or 3)
 */
class Diode
    :
    public Simget
{
    typedef Simget super;
public:

    //! Polynomial order \f$P_2\f$
    static const uint16_type Order = 1;
    static const uint16_type OrderGeo = 1;

    //! numerical type is double
    typedef double value_type;

    //! linear algebra backend factory
    typedef Backend<value_type> backend_type;
    //! linear algebra backend factory shared_ptr<> type
    typedef boost::shared_ptr<backend_type> backend_ptrtype;


    //! sparse matrix type associated with backend
    typedef backend_type::sparse_matrix_type sparse_matrix_type;
    //! sparse matrix type associated with backend (shared_ptr<> type)
    typedef backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    //! vector type associated with backend
    typedef backend_type::vector_type vector_type;
    //! vector type associated with backend (shared_ptr<> type)
    typedef backend_type::vector_ptrtype vector_ptrtype;

    //! geometry entities type composing the mesh, here Simplex in Dimension Dim of Order 2
    typedef Simplex<2,OrderGeo> convex_type;
    typedef Simplex<2,1> convex1_type;
    //typedef Hypercube<2,Order> convex_type;
    //! mesh type
    typedef Mesh<convex_type> mesh_type;
    //! mesh shared_ptr<> type
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    //! the basis type of our approximation space
    typedef bases<Lagrange<Order,Scalar,Discontinuous> > basis_type;


    //! the approximation function space type
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    //! the approximation function space type (shared_ptr<> type)
    typedef boost::shared_ptr<space_type> space_ptrtype;
    //! an element type of the approximation function space
    typedef space_type::element_type element_type;

    typedef bases<Lagrange<Order,Scalar> > c_basis_type;
    typedef FunctionSpace<mesh_type, c_basis_type> c_space_type;
    typedef boost::shared_ptr<c_space_type> c_space_ptrtype;
    typedef c_space_type::element_type c_element_type;

    //! the exporter factory type
    typedef Exporter<mesh_type,OrderGeo> export_type;
    //! the exporter factory (shared_ptr<> type)
    typedef boost::shared_ptr<export_type> export_ptrtype;

    /**
     * Constructor
     */
    Diode( po::variables_map const& vm, AboutData const& about )
        :
        super( vm, about ),
        M_backend( backend_type::build( this->vm() ) ),
        meshSize( this->vm()["hsize"].as<double>() ),
        timeStep( this->vm()["dt"].as<double>() ),
        Tfinal( this->vm()["Tfinal"].as<double>() ),
        rkmethod( (RKMethod)this->vm()["rkmethod"].as<int>() ),
        convex( this->vm()["convex"].as<std::string>() )
    {
    }

    template<typename ExExpr> void checkDG( ExExpr expr );

    void run();

    void run( const double* X, unsigned long P, double* Y, unsigned long N );

    template<typename BdyExpr>
    void
    FSolve( BdyExpr wbdy,
           element_type const& Ex, element_type const& Ey,element_type const& Bz,
           element_type& Exstar, element_type& Eystar, element_type& Bzstar );

    template<typename BdyExpr>
    void
    FSolve2( BdyExpr wbdy,
           element_type const& Ex, element_type const& Ey,element_type const& Bz,
           element_type& Exstar, element_type& Eystar, element_type& Bzstar );

    template<typename BdyExpr>
    void
    EulerStep(double& time, double dt, BdyExpr& wbdy,
                  element_type& Ex, element_type& Ey,element_type& Bz);
    template<typename BdyExpr>
    void
    EulerModifiedStep(double& time, double dt, BdyExpr& wbdy,
                  element_type& Ex, element_type& Ey,element_type& Bz);
    template<typename BdyExpr>
    void
    HeunStep(double& time, double dt, BdyExpr& wbdy,
                  element_type& Ex, element_type& Ey,element_type& Bz);
    template<typename BdyExpr>
    void
    RK4Step(double& time, double dt, BdyExpr& wbdy,
                  element_type& Ex, element_type& Ey,element_type& Bz);
private:

    //! linear algebra backend
    backend_ptrtype M_backend;

    //! mesh characteristic size
    double meshSize;
    double timeStep;
    double Tfinal;
    RKMethod rkmethod;
    //! convex of the domain
    std::string convex;

    sparse_matrix_ptrtype D;
    vector_ptrtype F_Ex, F_Ey, F_Bz;

    mesh_ptrtype mesh;
    space_ptrtype Xh;
}; // Diode
const uint16_type Diode::Order;
const uint16_type Diode::OrderGeo;

void
Diode::run()
{
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Execute Diode\n";
    std::vector<double> X( 2 );
    X[0] = meshSize;
    std::vector<double> Y( 3 );
    run( X.data(), X.size(), Y.data(), Y.size() );
}
template<typename ExExpr>
void
Diode::checkDG( ExExpr expr )
{
    // check continuity
    auto ijump  = integrate( internalfaces(mesh), trans(jumpv(expr))*jumpv(expr)  ).evaluate()( 0, 0 );
    std::cout << "continuity 1:" <<  math::sqrt( ijump ) << "\n";
#if 1
    auto v = vec(expr,expr,expr);
    auto vL = vec(leftfacev(expr),leftfacev(expr),leftfacev(expr));
    auto vR = vec(rightfacev(expr),rightfacev(expr),rightfacev(expr));
    auto A1 = trans(vec( cst(1.), cst(1.), cst(1.0) ));
    auto A2 = trans(vec( cst(1.), cst(1.), cst(1.0) ));
    auto myjump1 = (A1*(vL-vR))*leftfacev(N());
    auto myjump2 = (A2*(vL-vR))*leftfacev(N());
    auto ijump1 = integrate( internalfaces(mesh), trans(myjump1)*myjump1  ).evaluate()( 0, 0 );
    auto ijump2 = integrate( internalfaces(mesh), trans(myjump2)*myjump2  ).evaluate()( 0, 0 );
    std::cout << "continuity 1:" << math::sqrt(ijump1) << "\n";;
    std::cout << "continuity 2:" << math::sqrt(ijump2) << "\n";
#endif
}


template<typename BdyExpr>
void
Diode::FSolve2( BdyExpr wbdy,
              element_type const& Ex, element_type const& Ey,element_type const& Bz,
              element_type& dtEx, element_type& dtEy, element_type& dtBz )

{
    //Solve dtw*M = l(W)
    //l(W) = int (Ai di phi.w + bord)
    auto w = vec(idv(Ex),idv(Ey),idv(Bz));
    auto wR = vec(rightfacev(idv(Ex)),rightfacev(idv(Ey)),rightfacev(idv(Bz)));
    auto wL = vec(leftfacev(idv(Ex)),leftfacev(idv(Ey)),leftfacev(idv(Bz)));
    auto lEx=form1( _test=Xh, _vector=F_Ex, _init=true );
    auto lEy=form1( _test=Xh, _vector=F_Ey, _init=true );
    auto lBz=form1( _test=Xh, _vector=F_Bz, _init=true );

    // auto Anp_1 = vec( +Ny() * Ny() / 0.2e1, -Nx() * Ny() / 0.2e1, -Ny() / 0.2e1 );
    // auto Anp_2 = vec(-Nx() * Ny() / 0.2e1, Nx() * Nx() / 0.2e1, Nx() / 0.2e1 );
    // auto Anp_3 = vec( -Ny() / 0.2e1, Nx() / 0.2e1, cst(0.1e1 / 0.2e1) );

    // auto Anm_1 = vec(  -Ny() * Ny() / 0.2e1, Nx() * Ny() / 0.2e1, -Ny() / 0.2e1 );
    // auto Anm_2 = vec( Nx() * Ny() / 0.2e1, -Nx() * Nx() / 0.2e1, Nx() / 0.2e1);
    // auto Anm_3 = vec( -Ny() / 0.2e1, Nx() / 0.2e1, cst(-0.1e1 / 0.2e1) );

    //Aini+
    auto Anp_1 = vec(Ny() * Ny() /2.0, -Nx() * Ny() / 2.0, -Ny() / 2.0 );
    auto Anp_2 = vec(-Nx() * Ny() / 2.0, Nx() * Nx() / 2.0, Nx() / 2.0 );
    auto Anp_3 = vec( -Ny() / 2.0, Nx() / 2.0, cst(0.5) );

    //Aini-
    auto Anm_1 = vec(  -Ny() * Ny() / 2.0, Nx() * Ny() /2.0, -Ny() / 2.0 );
    auto Anm_2 = vec( Nx() * Ny() / 2.0, -Nx() * Nx() / 2.0, Nx() / 2.0);
    auto Anm_3 = vec( -Ny() / 2.0, Nx() /2.0, cst(0.5) );

    auto u = Xh->element();

    // update right hand side
    lEx = integrate( elements( mesh ), -dyv(Ex)*id(u));
    lEx += integrate( internalfaces(mesh),
                      trans(Anm_1)*(wR-wL)*rightface(id(u)) +
                      trans(Anp_1)*(wR-wL)*leftface(id(u)) );
    lEx += integrate( boundaryfaces(mesh), trans(Anp_1)*(wbdy-w)*id(u) );

    lEy = integrate( elements( mesh ), dxv(Ey)*id(u));
    lEy += integrate( internalfaces(mesh),
                      trans(Anm_2)*(wR-wL)*rightface(id(u)) +
                      trans(Anp_2)*(wR-wL)*leftface(id(u)) );
    lEy += integrate( boundaryfaces(mesh), trans(Anp_2)*(wbdy-w)*id(u) );

    lBz = integrate( elements( mesh ), (dxv(Bz) - dyv(Bz))*id(u));
    lBz += integrate( internalfaces(mesh),
                      trans(Anm_3)*(wR-wL)*rightface(id(u)) +
                      trans(Anp_3)*(wR-wL)*leftface(id(u)) );
    lBz += integrate( boundaryfaces(mesh), trans(Anp_3)*(wbdy-w)*id(u) );

    M_backend->solve( _matrix=D, _solution=dtEx, _rhs=F_Ex  );
    M_backend->solve( _matrix=D, _solution=dtEy, _rhs=F_Ey  );
    M_backend->solve( _matrix=D, _solution=dtBz, _rhs=F_Bz  );

}

template<typename BdyExpr>
void
Diode::FSolve( BdyExpr wbdy,
               element_type const& Ex, element_type const& Ey,element_type const& Bz,
               element_type& dtEx, element_type& dtEy, element_type& dtBz )

{
    //Solve dtw*M = l(W)
    //l(W) = int (Ai di phi.w + bord)
    auto w = vec(idv(Ex),idv(Ey),idv(Bz));
    auto wR = vec(rightfacev(idv(Ex)),rightfacev(idv(Ey)),rightfacev(idv(Bz)));
    auto wL = vec(leftfacev(idv(Ex)),leftfacev(idv(Ey)),leftfacev(idv(Bz)));
    auto lEx=form1( _test=Xh, _vector=F_Ex, _init=true );
    auto lEy=form1( _test=Xh, _vector=F_Ey, _init=true );
    auto lBz=form1( _test=Xh, _vector=F_Bz, _init=true );

    // auto Anp_1 = vec( +Ny() * Ny() / 0.2e1, -Nx() * Ny() / 0.2e1, -Ny() / 0.2e1 );
    // auto Anp_2 = vec(-Nx() * Ny() / 0.2e1, Nx() * Nx() / 0.2e1, Nx() / 0.2e1 );
    // auto Anp_3 = vec( -Ny() / 0.2e1, Nx() / 0.2e1, cst(0.1e1 / 0.2e1) );

    // auto Anm_1 = vec(  -Ny() * Ny() / 0.2e1, Nx() * Ny() / 0.2e1, -Ny() / 0.2e1 );
    // auto Anm_2 = vec( Nx() * Ny() / 0.2e1, -Nx() * Nx() / 0.2e1, Nx() / 0.2e1);
    // auto Anm_3 = vec( -Ny() / 0.2e1, Nx() / 0.2e1, cst(-0.1e1 / 0.2e1) );

    //Aini+
    auto Anp_1 = vec(Ny() * Ny() /2.0, -Nx() * Ny() / 2.0, -Ny() / 2.0 );
    auto Anp_2 = vec(-Nx() * Ny() / 2.0, Nx() * Nx() / 2.0, Nx() / 2.0 );
    auto Anp_3 = vec( -Ny() / 2.0, Nx() / 2.0, cst(0.5) );

    //Aini-
    auto Anm_1 = vec(  -Ny() * Ny() / 2.0, Nx() * Ny() /2.0, -Ny() / 2.0 );
    auto Anm_2 = vec( Nx() * Ny() / 2.0, -Nx() * Nx() / 2.0, Nx() / 2.0);
    auto Anm_3 = vec( -Ny() / 2.0, Nx() /2.0, cst(0.5) );

    auto u = Xh->element();

    // update right hand side
    lEx = integrate( elements( mesh ), -idv(Ex)*dy(u));
    std::cout<< "jump term Ex : "<<integrate(internalfaces(mesh),
                                             trans(Anm_1)*(wR-wL)*rightfacev(idv(Ex)) +
                                             trans(Anp_1)*(wR-wL)*leftfacev(idv(Ex))).evaluate()(0,0)<< std::endl;
    lEx += integrate( internalfaces(mesh),
                      trans(Anm_1)*(wR-wL)*rightface(id(u)) +
                      trans(Anp_1)*(wR-wL)*leftface(id(u)) );
    lEx += integrate( boundaryfaces(mesh), trans(Anp_1)*(wbdy-w)*id(u) );


    auto F_Ex1 = M_backend->newVector( Xh );
    auto lEx1=form1( _test=Xh, _vector=F_Ex1, _init=true );
    lEx1 = integrate( internalfaces(mesh),
                      trans(Anm_1)*(wR-wL)*rightface(id(u)) +
                      trans(Anp_1)*(wR-wL)*leftface(id(u)) );

    F_Ex1->close();
    F_Ex1->printMatlab( "jump_Ex.m" );
    std::cout << "internalfaces( Ex ) = " << Feel::dot( *F_Ex1, Ex ) << "\n";

    std::cout << "boundary term Ex : " <<integrate( boundaryfaces(mesh), trans(Anm_1)*(wbdy-w)*idv(Ex) ).evaluate() << endl;

    lEx1 = integrate( boundaryfaces(mesh), trans(Anm_1)*(wbdy-w)*id(u) );
    F_Ex1->close();
    F_Ex1->printMatlab( "bdy_Ex.m" );
    std::cout << "internalfaces( Ex ) = " << Feel::dot( *F_Ex1, Ex ) << "\n";

    // Ey
    lEy = integrate( elements( mesh ), idv(Ey)*dx(u));

    std::cout<< "jump term Ey : "<<integrate(internalfaces(mesh),
                                             trans(Anm_2)*(wR-wL)*rightfacev(idv(Ey)) +
                                             trans(Anp_2)*(wR-wL)*leftfacev(idv(Ey))).evaluate()(0,0)<< std::endl;
    lEy += integrate( internalfaces(mesh),
                      trans(Anm_2)*(wR-wL)*rightface(id(u)) +
                      trans(Anp_2)*(wR-wL)*leftface(id(u)) );
    lEy += integrate( boundaryfaces(mesh), trans(Anp_2)*(wbdy-w)*id(u) );

    lEx1 = integrate( internalfaces(mesh),
                      trans(Anm_2)*(wR-wL)*rightfacev(idv(Ey)) +
                      trans(Anp_2)*(wR-wL)*leftfacev(idv(Ey)));
    F_Ex1->close();
    std::cout << "internalfaces( Ex ) = " << Feel::dot( *F_Ex1, Ey ) << "\n";
    std::cout << "boundary term Ey : " <<integrate( boundaryfaces(mesh), trans(Anm_2)*(wbdy-w)*idv(Ey) ).evaluate() << endl;


    lBz = integrate( elements( mesh ), idv(Bz)*(dx(u) - dy(u)));


    std::cout<< "jump term Bz : "<<integrate(internalfaces(mesh),
                                             trans(Anm_3)*(wR-wL)*rightfacev(idv(Bz)) +
                                             trans(Anp_3)*(wR-wL)*leftfacev(idv(Bz))).evaluate()(0,0)<< std::endl;
    lBz += integrate( internalfaces(mesh),
                      trans(Anm_3)*(wR-wL)*rightface(id(u)) +
                      trans(Anp_3)*(wR-wL)*leftface(id(u)) );
    lBz += integrate( boundaryfaces(mesh), trans(Anp_3)*(wbdy-w)*id(u) );

    std::cout << "boundary term Bz : " <<integrate( boundaryfaces(mesh), trans(Anm_3)*(wbdy-w)*idv(Bz) ).evaluate() << endl;
    lEx1 = integrate( boundaryfaces(mesh), trans(Anm_3)*(wbdy-w)*id(u) );
    F_Ex1->close();
    F_Ex1->printMatlab( "bdy_Bz.m" );

    std::cout << "bdy( Bz ) = " << Feel::dot( *F_Ex1, Bz ) << "\n";

    F_Ex->printMatlab( "F_Ex.m" );
    F_Ey->printMatlab( "F_Ey.m" );
    F_Bz->printMatlab( "F_Bz.m" );
    M_backend->solve( _matrix=D, _solution=dtEx, _rhs=F_Ex  );
    M_backend->solve( _matrix=D, _solution=dtEy, _rhs=F_Ey  );
    M_backend->solve( _matrix=D, _solution=dtBz, _rhs=F_Bz  );

}
template<typename BdyExpr>
void
Diode::EulerStep(double& time,double dt,
                 BdyExpr& wbdy,
                 element_type& Ex, element_type& Ey,element_type& Bz )
{
    element_type Exstar = Xh->element();
    element_type Eystar = Xh->element();
    element_type Bzstar = Xh->element();
    FSolve(wbdy, Ex, Ey, Bz, Exstar, Eystar, Bzstar);

    Ex.add(dt, Exstar);
    Ey.add(dt, Eystar);
    Bz.add(dt, Bzstar);
    time += dt;
}
template<typename BdyExpr>
void
Diode::EulerModifiedStep(double& time,double dt,
                         BdyExpr& wbdy,
                         element_type& Ex, element_type& Ey,element_type& Bz )
{
    element_type Exstar = Xh->element();
    element_type Eystar = Xh->element();
    element_type Bzstar = Xh->element();
    element_type Exn = Ex;
    element_type Eyn = Ey;
    element_type Bzn = Bz;
    FSolve(wbdy, Ex, Ey, Bz, Exstar, Eystar, Bzstar);

    Exn.add(dt/2.0, Exstar);
    Eyn.add(dt/2.0, Eystar);
    Bzn.add(dt/2.0, Bzstar);
    time += dt/2.0;

    FSolve(wbdy, Exn, Eyn, Bzn, Exstar, Eystar, Bzstar);

    Ex.add(dt, Exstar);
    Ey.add(dt, Eystar);
    Bz.add(dt, Bzstar);
    time += dt/2.0;

}
template<typename BdyExpr>
void
Diode::HeunStep(double& time,double dt,
                         BdyExpr& wbdy,
                         element_type& Ex, element_type& Ey,element_type& Bz )
{
    element_type Exstar = Xh->element();
    element_type Eystar = Xh->element();
    element_type Bzstar = Xh->element();
    element_type Exn = Ex;
    element_type Eyn = Ey;
    element_type Bzn = Bz;
    FSolve(wbdy, Ex, Ey, Bz, Exstar, Eystar, Bzstar);

    Ex.add(dt/2.0, Exstar);
    Ey.add(dt/2.0, Eystar);
    Bz.add(dt/2.0, Bzstar);
    Exn.add(dt, Exstar);
    Eyn.add(dt, Eystar);
    Bzn.add(dt, Bzstar);
    time += dt;

    FSolve(wbdy, Exn, Eyn, Bzn, Exstar, Eystar, Bzstar);

    Ex.add(dt/2.0, Exstar);
    Ey.add(dt/2.0, Eystar);
    Bz.add(dt/2.0, Bzstar);

}
template<typename BdyExpr>
void
Diode::RK4Step(double& time,double dt,
                         BdyExpr& wbdy,
                         element_type& Ex, element_type& Ey,element_type& Bz )
{
    element_type Exstar = Xh->element();
    element_type Eystar = Xh->element();
    element_type Bzstar = Xh->element();
    element_type Exn = Ex;
    element_type Eyn = Ey;
    element_type Bzn = Bz;
    element_type Exk = Ex;
    element_type Eyk = Ey;
    element_type Bzk = Bz;

    //k1
    FSolve(wbdy, Ex, Ey, Bz, Exk, Eyk, Bzk);

    Ex.add(dt/6.0, Exk);
    Ey.add(dt/6.0, Eyk);
    Bz.add(dt/6.0, Bzk);
    Exstar=Exn;
    Eystar=Eyn;
    Bzstar=Bzn;
    Exstar.add(dt/2.0, Exk);
    Eystar.add(dt/2.0, Eyk);
    Bzstar.add(dt/2.0, Bzk);

    //k2
    time += dt/2.0;
    FSolve(wbdy, Exstar, Eystar, Bzstar, Exk, Eyk, Bzk);

    Ex.add(dt/3.0, Exk);
    Ey.add(dt/3.0, Eyk);
    Bz.add(dt/3.0, Bzk);
    Exstar=Exn;
    Eystar=Eyn;
    Bzstar=Bzn;
    Exstar.add(dt/2.0, Exk);
    Eystar.add(dt/2.0, Eyk);
    Bzstar.add(dt/2.0, Bzk);

    //k3
    FSolve(wbdy, Exstar, Eystar, Bzstar, Exk, Eyk, Bzk);

    Ex.add(dt/3.0, Exk);
    Ey.add(dt/3.0, Eyk);
    Bz.add(dt/3.0, Bzk);
    Exstar=Exn;
    Eystar=Eyn;
    Bzstar=Bzn;
    Exstar.add(dt, Exk);
    Eystar.add(dt, Eyk);
    Bzstar.add(dt, Bzk);

    //k4
    time += dt/2.0;
    FSolve(wbdy, Exstar, Eystar, Bzstar, Exk, Eyk, Bzk);

    Ex.add(dt/6.0, Exk);
    Ey.add(dt/6.0, Eyk);
    Bz.add(dt/6.0, Bzk);

}

void
Diode::run( const double* X, unsigned long P, double* Y, unsigned long N )
{
    if ( !this->vm().count( "nochdir" ) )
        Environment::changeRepository( boost::format( "examples/maxwell/%1%/%2%/P%3%/h_%4%/" )
                                       % this->about().appName()
                                       % convex
                                       % OrderGeo
                                       % X[0] );
    mesh = createGMSHMesh( _mesh=new mesh_type,
                           _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                           _desc=diodegeo(X[0],OrderGeo,convex) );

    /**
     * The function space and some associated elements(functions) are then defined
     */
    /** \code */
    Xh = space_type::New( mesh );
    auto Xhc = c_space_type::New( mesh );
    auto Ex = Xh->element();
    auto Ey = Xh->element();
    auto Bz = Xh->element();
    auto dtEx = Xh->element();
    auto dtEy = Xh->element();
    auto dtBz = Xh->element();
    auto dtEx2 = Xh->element();
    auto dtEy2 = Xh->element();
    auto dtBz2 = Xh->element();
    auto Exe = Xhc->element();
    auto Eye = Xhc->element();
    auto Bze = Xhc->element();
    auto Exc = Xhc->element();
    auto Eyc = Xhc->element();
    auto Bzc = Xhc->element();
    auto u = Xh->element();
    auto v = Xh->element();

    using namespace Feel::vf;
    double dt=std::min(0.1*meshSize/(Order+1),timeStep);

    double pi=M_PI;
    double k=2*pi; //=0;
    double theta=pi/4;
    //theta=0;
    double vu=cos(theta);
    double vv=sin(theta);
    double time = 0;
    auto c=cos(k * (vu * Px() + vv * Py() - cst_ref(time)));
    auto Ex_exact = -vv*c;
    auto Ey_exact = vu*c;
    auto Bz_exact = c;
    auto w_exact = vec(Ex_exact, Ey_exact, Bz_exact );
    auto wL_exact = vec(leftfacev(Ex_exact), leftfacev(Ey_exact), leftfacev(Bz_exact) );
    auto wR_exact = vec(rightfacev(Ex_exact), rightfacev(Ey_exact), rightfacev(Bz_exact) );

    checkDG( Ex_exact );
    checkDG( Ey_exact );
    checkDG( Bz_exact );

    F_Ex = M_backend->newVector( Xh );
    F_Ey = M_backend->newVector( Xh );
    F_Bz = M_backend->newVector( Xh );

    // left hand side
    D=M_backend->newMatrix( Xh, Xh );
    auto a = form2( _test=Xh, _trial=Xh, _matrix=D, _init=true );
    a = integrate( elements(mesh), idt(Ex)*id(u) );
    D->printMatlab("mass.m");
    auto backend = backend_type::build( this->vm() );
    auto exporter( export_type::New( this->vm(),
                                     (boost::format( "%1%-%2%" )
                                      % this->about().appName()
                                      % convex).str() ) );
    Ex = vf::project( Xh, elements(mesh), Ex_exact );
    Ey = vf::project( Xh, elements(mesh), Ey_exact );
    Bz = vf::project( Xh, elements(mesh), Bz_exact );

    auto L2Proj = projector( Xhc, Xhc );
    Exc = L2Proj->project( idv(Ex) );
    Eyc = L2Proj->project( idv(Ey) );
    Bzc = L2Proj->project( idv(Bz) );

    Exe = vf::project( Xhc, elements(mesh), Ex_exact );
    Eye = vf::project( Xhc, elements(mesh), Ey_exact );
    Bze = vf::project( Xhc, elements(mesh), Bz_exact );

    exporter->step(time)->setMesh( mesh );
#if 0
    exporter->step(time)->add( "Ex", Ex );
    exporter->step(time)->add( "Ey", Ey );
    exporter->step(time)->add( "Bz", Bz );
#endif
    exporter->step(time)->add( "Exc", Exc );
    exporter->step(time)->add( "Eyc", Eyc );
    exporter->step(time)->add( "Bzc", Bzc );
    exporter->step(time)->add( "ExExact", Exe );
    exporter->step(time)->add( "EyExact", Eye );
    exporter->step(time)->add( "BzExact", Bze );

    auto w = vec(idv(Ex),idv(Ey),idv(Bz));
    auto wR = vec(rightfacev(idv(Ex)),rightfacev(idv(Ey)),rightfacev(idv(Bz)));
    auto wL = vec(leftfacev(idv(Ex)),leftfacev(idv(Ey)),leftfacev(idv(Bz)));

    std::cout<<rkmethod<<endl;
    switch( rkmethod )
        {
        case EULER_EXPLICIT:
            std::cout << "Euler explicit" << endl;
            break;
        case EULER_MODIFIED:
            std::cout << "Euler modified" << endl;
            break;
        case HEUN:
            std::cout << "Heun" << endl;
            break;
        case RK4:
            std::cout << "RK4" << endl;
            break;
        }
    exporter->save();
    std::cout << "Saved initial/exact solution\n";
    for( time = 0; time <= Tfinal; )
    {
        std::cout << "============================================================" << std::endl;
        std::cout << "time = " << time << "s, dt=" << dt << ", final time=" << Tfinal << std::endl;

        Ex = vf::project( Xh, elements(mesh), Ex_exact );
        Ey = vf::project( Xh, elements(mesh), Ey_exact );
        Bz = vf::project( Xh, elements(mesh), Bz_exact );

        checkDG( idv(Ex) );
        checkDG( idv(Ey) );
        checkDG( idv(Bz) );
        FSolve(w, Ex, Ey, Bz, dtEx, dtEy, dtBz);
        FSolve2(w, Ex, Ey, Bz, dtEx2, dtEy2, dtBz2);
        time += dt;

        std::cout << "||exact-Ex||_2" << integrate(elements(mesh), (idv(dtEx)+(idv(Ex)-Ex_exact)/dt)*(idv(dtEx)+(idv(Ex)-Ex_exact)/dt) ).evaluate().norm() << std::endl;
        std::cout << "||exact-Ey||_2" << integrate(elements(mesh), (idv(dtEy)+(idv(Ey)-Ey_exact)/dt)*(idv(dtEy)+(idv(Ey)-Ey_exact)/dt) ).evaluate().norm() << std::endl;
        std::cout << "||exact-Bz||_2" << integrate(elements(mesh), (idv(dtBz)+(idv(Bz)-Bz_exact)/dt)*(idv(dtBz)+(idv(Bz)-Bz_exact)/dt) ).evaluate().norm() << std::endl;

        std::cout << "||dtEx-dtEx2||_2" << integrate(elements(mesh), (idv(dtEx)-idv(dtEx2))*(idv(dtEx)-idv(dtEx2) )).evaluate().norm() << std::endl;
        std::cout << "||dtEy-dtEy2||_2" << integrate(elements(mesh), (idv(dtEy)-idv(dtEy2))*(idv(dtEy)-idv(dtEy2) )).evaluate().norm() << std::endl;
        std::cout << "||dtBz-dtBz2||_2" << integrate(elements(mesh), (idv(dtBz)-idv(dtBz2))*(idv(dtBz)-idv(dtBz2) )).evaluate().norm() << std::endl;
        /*
        switch( rkmethod )
            {
            case EULER_EXPLICIT:
                EulerStep(time, dt, w_exact, Ex, Ey, Bz);
                break;
            case EULER_MODIFIED:
                EulerModifiedStep(time, dt, w_exact, Ex, Ey, Bz);
                break;
            case HEUN:
                HeunStep(time, dt, w_exact, Ex, Ey, Bz);
                break;
            case RK4:
                RK4Step(time, dt, w_exact, Ex, Ey, Bz);
                break;
            }


        std::cout << "||exact-Ex||_2 = " << integrate(elements(mesh), (idv(Ex)-Ex_exact)*(idv(Ex)-Ex_exact) ).evaluate().norm() << std::endl;
        std::cout << "||exact-Ey||_2 = " << integrate(elements(mesh), (idv(Ey)-Ey_exact)*(idv(Ey)-Ey_exact) ).evaluate().norm() << std::endl;
        std::cout << "||exact-Bz||_2 = " << integrate(elements(mesh), (idv(Bz)-Bz_exact)*(idv(Bz)-Bz_exact) ).evaluate().norm() << std::endl;*/

        // save
        exporter->step(time)->setMesh( mesh );
#if 0
        exporter->step(time)->add( "Ex", Ex );
        exporter->step(time)->add( "Ey", Ey );
        exporter->step(time)->add( "Bz", Bz );
#endif
        Exc = L2Proj->project( idv(Ex) );
        Eyc = L2Proj->project( idv(Ey) );
        Bzc = L2Proj->project( idv(Bz) );
        exporter->step(time)->add( "Exc", Exc );
        exporter->step(time)->add( "Eyc", Eyc );
        exporter->step(time)->add( "Bzc", Bzc );

        Exe = vf::project( Xhc, elements(mesh), Ex_exact );
        Eye = vf::project( Xhc, elements(mesh), Ey_exact );
        Bze = vf::project( Xhc, elements(mesh), Bz_exact );
        exporter->step(time)->add( "ExExact", Exe );
        exporter->step(time)->add( "EyExact", Eye );
        exporter->step(time)->add( "BzExact", Bze );

        exporter->save();
    }



} // Diode::run

/**
 * main function: entry point of the program
 */
int
main( int argc, char** argv )
{
    Environment env( argc, argv );
    /**
     * create an application
     */
    /** \code */
    Application app( argc, argv, makeAbout(), makeOptions() );
    if ( app.vm().count( "help" ) )
    {
        std::cout << app.optionsDescription() << "\n";
        return 0;
    }
    /** \endcode */

    /**
     * register the simgets
     */
    /** \code */
    app.add( new Diode( app.vm(), app.about() ) );
    /** \endcode */

    /**
     * run the application
     */
    /** \code */
    app.run();
    /** \endcode */
}






