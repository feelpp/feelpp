/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
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
 \page Diode 2D Maxwell Dimulation in a Diode
\author Thomas Strub
\author Philippe Helluy
\author Christophe Prud'homme
\date 2011-06-01

\tableofcontents
<br>
<hr>
<br>

\section Diode_Description Description
The Maxwell equations read:
\f{eqnarray*}
\frac{-1}{c^{2}}\frac{\partial {{\bm E}}}{\partial t}+\nabla\times {{\bm B}}& = & \mu_{0} {{\bm J}}\\
{{\bm B}}{t}+\nabla\times {{\bm E}}& = & 0\\
\nabla \cdot {{\bm B}}& = & 0\\
\nabla \cdot {{\bm E}}& = & \frac{\rho}{\epsilon_{o}}
\f}

<!--
\f{eqnarray*}
\frac{-1}{c^{2}}\frac{\partial \ensuremath{{\bm E}}\xspace}{\partial t}+\nabla\times \ensuremath{{\bm B}}\xspace & = & \mu_{0} \ensuremath{{\bm J}}\xspace\\
\ensuremath{{\bm B}}\xspace_{t}+\nabla\times \ensuremath{{\bm E}}\xspace & = & 0\\
\nabla \cdot \ensuremath{{\bm B}}\xspace & = & 0\\
\nabla \cdot \ensuremath{{\bm E}}\xspace & = & \frac{\rho}{\epsilon_{o}}
\f}
-->

<!--
where \f$\ensuremath{{\bm E}}\xspace\f$ is the electric field, \f$\ensuremath{{\bm B}}\xspace\f$ the magnetic field,
\f$\ensuremath{{\bm J}}\xspace\f$ the current density, \f$ c \f$ the speed of light, \f$ \rho \f$
-->
where \f${{\bm E}}\f$ is the electric field, \f${{\bm B}}\f$ the magnetic field,
\f${{\bm J}}\f$ the current density, \f$ c \f$ the speed of light, \f$ \rho \f$
density of electric charge, \f$ \mu_ {0} \f$ the vacuum permeability
and \f$ \epsilon_ {0} \f$ the vacuum permittivity.<br>
In the midst industrial notament in aeronautics, systems
Products must verify certain standards such as the receipt
an electromagnetic wave emitted by a radar does not cause
the inefficassité of part or all of the hardware in the
system.<br>
Thus, the simulation of such situations can develop when
or during the certification of a new product to test its reaction
to such attacks.<br>
Also note that the last two equations are actually
initial conditions, since if we assume they are true at the moment
\f$ t = 0\f$ then it can be deduced from the first two.<br>
At \f$t=0s\f$, we suppose that
\f{eqnarray}
\nabla \cdot {{\bm B}}& = & 0\\
\nabla \cdot {{\bm E}}& = & \frac{\rho}{\epsilon_{o}}
\f}
<!--
\f{eqnarray}
\nabla \cdot \ensuremath{{\bm B}}\xspace & = & 0\\
\nabla \cdot \ensuremath{{\bm E}}\xspace & = & \frac{\rho}{\epsilon_{o}}
\f}
-->

Suppose that \f${{\bm B}} = (B_x, B_y, B_z )^T\f$ and \f${{\bm E}}=(E_x,E_y,E_z)^T\f$
<!-- Suppose that \f$\ensuremath{{\bm B}}\xspace = (B_x, B_y, B_z )^T\f$ and \f$\ensuremath{{\bm E}}\xspace=(E_x,E_y,E_z)^T\f$ -->
i.e.\f{eqnarray}
\frac{\partial B_{x}}{\partial x}(t=0)+\frac{\partial B_{y}}{\partial y}(t=0)+\frac{\partial B_{z}}{\partial z} & (t=0)= & 0\\
\frac{\partial E_{x}}{\partial x}(t=0)+\frac{\partial E_{y}}{\partial
y}(t=0)+\frac{\partial E_{z}}{\partial z} & (t=0)= & \frac{\rho}{\epsilon_{o}}
\f}

Differentiating the first of these two equations with respect to time,
we get:
\f{multline}
  \label{eq:6}
\frac{\partial}{\partial t}\frac{\partial}{\partial
x}B_{x}+\frac{\partial}{\partial t}\frac{\partial}{\partial y}B_{y}+\frac{\partial}{\partial t}\frac{\partial}{\partial z}B_{z} = \\
\frac{\partial}{\partial x}\left(\frac{\partial}{\partial y}E_{z}-\frac{\partial}{\partial z}E_{y}\right)+\frac{\partial}{\partial y}\left(\frac{\partial}{\partial z}E_{x}-\frac{\partial}{\partial x}E_{z}\right)+\frac{\partial}{\partial z}\left(\frac{\partial}{\partial x}E_{y}-\frac{\partial}{\partial y}E_{x}\right)\\
 =  0
\f}

thanks to
\f{equation}
  \label{eq:3}
  {{\bm B}}_{t}+\nabla\times {{\bm E}} =0
\f}
<!--
\f{equation}
  \label{eq:3}
  \ensuremath{{\bm B}}\xspace_{t}+\nabla\times \ensuremath{{\bm E}}\xspace =0
\f}
-->

So, for all \f$t\geq0\f$,
\f{equation}
  \label{eq:4}
  \nabla \cdot {{\bm B}}(t)=\nabla\cdot(0)=0
\f}

<!--
So, for all \f$t\geq0\f$,
\f{equation}
  \label{eq:4}
  \nabla \cdot \ensuremath{{\bm B}}\xspace(t)=\nabla\cdot\B(0)=0
\f}
-->

We deduce the same way the second equation, using the charge
conservation equation :
\f{equation}
  \label{eq:2}
  \frac{\partial \rho}{\partial t} + \nabla \cdot  (\rho {{\bm J}}) = 0
\f}
<!--
\f{equation}
  \label{eq:2}
  \frac{\partial \rho}{\partial t} + \nabla \cdot  (\rho \ensuremath{{\bm J}}\xspace) = 0
\f}
-->

\section Diode_Theory Theory


\section Diode_Implementation Implementation


\section Diode_Results Results


*/

#include <feel/options.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
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
    po::options_description diodeoptions( "Diode options" );
    diodeoptions.add_options()
    ( "verbose", po::value<bool>()->default_value( 1 ), "verbose output" )
    ( "geomap", po::value<int>()->default_value( 0 ), "type of geomap" )
    ( "hsize", po::value<double>()->default_value( 0.5 ), "mesh size" )
    ( "convex", Feel::po::value<std::string>()->default_value( "simplex" ), "shape of the convex used in the mesh (either simplex or hypercube)" )
    ( "penaldir", Feel::po::value<double>()->default_value( 20 ), "penalisation parameter for the weak boundary conditions" )
    ( "dt", Feel::po::value<double>()->default_value( 0.01 ), "timestep value" )
    ( "Tfinal", Feel::po::value<double>()->default_value( 1 ), "final time" )
    ( "rkmethod", Feel::po::value<int>()->default_value( 2 ), "rk method, 0=euler, 1=modified euler, 2=heun, 3=rk4" )
    ( "theta", Feel::po::value<double>()->default_value( 0.0 ), "angle of propagation" )
    ( "metalFlag", Feel::po::value<int>()->default_value( 1 ), "specify the flag of faces with metalic condition" )
    ( "fieldFlag", Feel::po::value<int>()->default_value( 2 ), "specify the flag of faces where apply a field" )
    ( "initfields", Feel::po::value<int>()->default_value( 1 ), "initialise the fields" )
    ( "exportfreq", Feel::po::value<int>()->default_value( 1 ), "results exporting frequency" )
    ;
    return diodeoptions.add( Feel::feel_options() ).add( backend_options( "mass" ) );
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
                     "Copyright (c) 2011 Universite Joseph Fourier" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Philippe Helluy", "developer", "helluy@math.unistra.fr", "" );
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
    static const uint16_type Order = 4;
    static const uint16_type OrderGeo = 4;

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
    Diode()
        :
        super(),
        M_backend( backend_type::build( "mass" ) ),
        verbose(  boption("verbose") ),
        meshSize( doption("hsize") ),
        timeStep( doption("dt") ),
        Tfinal(   doption("Tfinal") ),
        rkmethod( ( RKMethod )ioption("rkmethod") ),
        convex( soption("convex") )
    {
    }

    template<typename ExExpr> void checkDG( ExExpr expr );

    void run();

    void run( const double* X, unsigned long P, double* Y, unsigned long N );

    template<typename BdyExpr>
    void
    FSolve( BdyExpr& wbdy, double dt,
            element_type const& Ex, element_type const& Ey,element_type const& Bz,
            element_type& Exstar, element_type& Eystar, element_type& Bzstar );

    template<typename BdyExpr>
    void
    EulerStep( double& time, double dt, BdyExpr& wbdy,
               element_type& Ex, element_type& Ey,element_type& Bz );
    template<typename BdyExpr>
    void
    EulerModifiedStep( double& time, double dt, BdyExpr& wbdy,
                       element_type& Ex, element_type& Ey,element_type& Bz );
    template<typename BdyExpr>
    void
    HeunStep( double& time, double dt, BdyExpr& wbdy,
              element_type& Ex, element_type& Ey,element_type& Bz );
    template<typename BdyExpr>
    void
    RK4Step( double& time, double dt, BdyExpr& wbdy,
             element_type& Ex, element_type& Ey,element_type& Bz );
private:

    //! linear algebra backend
    backend_ptrtype M_backend;

    bool verbose;

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
    auto ijump  = integrate( internalfaces( mesh ), trans( jumpv( expr ) )*jumpv( expr )  ).evaluate()( 0, 0 );
    std::cout << "continuity 1:" <<  math::sqrt( ijump ) << "\n";
#if 1
    auto v = vec( expr,expr,expr );
    auto vL = vec( leftfacev( expr ),leftfacev( expr ),leftfacev( expr ) );
    auto vR = vec( rightfacev( expr ),rightfacev( expr ),rightfacev( expr ) );
    auto A1 = trans( vec( cst( 1. ), cst( 1. ), cst( 1.0 ) ) );
    auto A2 = trans( vec( cst( 1. ), cst( 1. ), cst( 1.0 ) ) );
    auto myjump1 = ( A1*( vL-vR ) )*leftfacev( N() );
    auto myjump2 = ( A2*( vL-vR ) )*leftfacev( N() );
    auto ijump1 = integrate( internalfaces( mesh ), trans( myjump1 )*myjump1  ).evaluate()( 0, 0 );
    auto ijump2 = integrate( internalfaces( mesh ), trans( myjump2 )*myjump2  ).evaluate()( 0, 0 );
    std::cout << "continuity 1:" << math::sqrt( ijump1 ) << "\n";;
    std::cout << "continuity 2:" << math::sqrt( ijump2 ) << "\n";
#endif
}


template<typename BdyExpr>
void
Diode::FSolve( BdyExpr& wbdy,
               double dt,
               element_type const& Exn,
               element_type const& Eyn,
               element_type const& Bzn,
               element_type& dtEx,
               element_type& dtEy,
               element_type& dtBz )

{
    boost::timer ti;
    //Solve dtw*M = l(W)
    //l(W) = int (Ai di phi.w + bord)
    auto w = vec( idv( Exn ),idv( Eyn ),idv( Bzn ) );
    auto wR = vec( rightfacev( idv( Exn ) ),rightfacev( idv( Eyn ) ),rightfacev( idv( Bzn ) ) );
    auto wL = vec( leftfacev( idv( Exn ) ),leftfacev( idv( Eyn ) ),leftfacev( idv( Bzn ) ) );
    auto wMetal = vec( -leftfacev( idv( Exn ) ),-leftfacev( idv( Eyn ) ),leftfacev( idv( Bzn ) ) );
    auto lEx=form1( _test=Xh, _vector=F_Ex, _init=true );
    auto lEy=form1( _test=Xh, _vector=F_Ey, _init=true );
    auto lBz=form1( _test=Xh, _vector=F_Bz, _init=true );

    auto Anp_1 = vec( +Ny() * Ny() / 0.2e1, -Nx() * Ny() / 0.2e1, -Ny() / 0.2e1 );
    auto Anp_2 = vec( -Nx() * Ny() / 0.2e1, Nx() * Nx() / 0.2e1, Nx() / 0.2e1 );
    auto Anp_3 = vec( -Ny() / 0.2e1, Nx() / 0.2e1, cst( 0.1e1 / 0.2e1 ) );

    auto Anm_1 = vec(  -Ny() * Ny() / 0.2e1, Nx() * Ny() / 0.2e1, -Ny() / 0.2e1 );
    auto Anm_2 = vec( Nx() * Ny() / 0.2e1, -Nx() * Nx() / 0.2e1, Nx() / 0.2e1 );
    auto Anm_3 = vec( -Ny() / 0.2e1, Nx() / 0.2e1, cst( -0.1e1 / 0.2e1 ) );

    auto v1 = vec( cst( 1.0 ) , cst( 0.0 ) , cst( 0.0 ) );
    auto v2 = vec( cst( 0.0 ) , cst( 1.0 ) , cst( 0.0 ) );
    auto v3 = vec( cst( 0.0 ) , cst( 0.0 ) , cst( 1.0 ) );

    auto u = Xh->element();

    //
    // Ex
    //
    lEx = integrate( _range=elements( mesh ),  _expr=id( u )*dyv( Bzn ) );
    lEx += integrate( _range=internalfaces( mesh ),
                      _expr=( trans( Anm_1 )*( wL-wR ) )*leftface( id( u ) )
                            + ( trans( Anp_1 )*( wL-wR ) )*rightface( id( u ) ) );
    //lEx += integrate( boundaryfaces(mesh), trans(Anm_1)*(wL-wbdy)*id(u) );
    // lEx += integrate( markedfaces(mesh, this->vm()["metalFlag"].as<int>() ), trans(Anm_1)*(wL-wMetal)*id(u));
    // lEx += integrate( markedfaces(mesh, this->vm()["fieldFlag"].as<int>() ), trans(v1)*wbdy*id(u));

    lEx += integrate( _range=markedfaces( mesh, "Metal" ), _expr=( trans( Anm_1 )*( wL-wMetal ) )*id( u ) );
    lEx += integrate( _range=markedfaces( mesh, "Dirichlet" ), _expr=( trans( Anm_1 )*( wL-wbdy ) )*id( u ) );

    //
    // Ey
    //
    lEy = integrate( _range=elements( mesh ),  _expr=-id( u )*dxv( Bzn ) );
    lEy += integrate( _range=internalfaces( mesh ),
                      _expr=( trans( Anm_2 )*( wL-wR ) )*leftface( id( u ) )
                            + ( trans( Anp_2 )*( wL-wR ) )*rightface( id( u ) ) );
    //lEy += integrate( boundaryfaces(mesh), trans(Anm_2)*(wL-wbdy)*id(u) );
    // lEy += integrate( markedfaces(mesh, this->vm()["metalFlag"].as<int>() ), trans(Anm_3)*(wL-wMetal)*id(u));
    // lEy += integrate( markedfaces(mesh, this->vm()["fieldFlag"].as<int>() ), trans(v2)*wbdy*id(u));

    lEy += integrate( _range=markedfaces( mesh, "Metal" ), _expr=( trans( Anm_2 )*( wL-wMetal ) )*id( u ) );
    lEy += integrate( _range=markedfaces( mesh, "Dirichlet" ), _expr=( trans( Anm_2 )*( wL-wbdy ) )*id( u ) );

    //
    // Bz
    //
    lBz = integrate( _range=elements( mesh ),  _expr=-id( u )*dxv( Eyn ) + id( u )*dyv( Exn ) );
    lBz += integrate( _range=internalfaces( mesh ),
                      _expr=( trans( Anm_3 )*( wL-wR ) )*leftface( id( u ) )
                            + ( trans( Anp_3 )*( wL-wR ) )*rightface( id( u ) ) );
    //lBz += integrate( boundaryfaces(mesh), trans(Anm_3)*(wL-wbdy)*id(u) );
    // lBz += integrate( markedfaces(mesh, this->vm()["metalFlag"].as<int>() ), trans(Anm_3)*(wL-wMetal)*id(u));
    // lBz += integrate( markedfaces(mesh, this->vm()["fieldFlag"].as<int>() ), trans(v3)*wbdy*id(u));

    lBz += integrate( _range=markedfaces( mesh, "Metal" ), _expr=( trans( Anm_3 )*( wL-wMetal ) )*id( u ) );
    lBz += integrate( _range=markedfaces( mesh, "Dirichlet" ), _expr=( trans( Anm_3 )*( wL-wbdy ) )*id( u ) );

    if ( verbose )
        std::cout << " -- assembly in " << ti.elapsed() << "s\n";

    ti.restart();

    //
    // Solve
    //
    M_backend->solve( _matrix=D, _solution=dtEx, _rhs=F_Ex  );
    M_backend->solve( _matrix=D, _solution=dtEy, _rhs=F_Ey  );
    M_backend->solve( _matrix=D, _solution=dtBz, _rhs=F_Bz  );

    if ( verbose )
        std::cout << " -- solve in " << ti.elapsed() << "s\n";
}

template<typename BdyExpr>
void
Diode::EulerStep( double& time,double dt,
                  BdyExpr& wbdy,
                  element_type& Ex, element_type& Ey,element_type& Bz )
{
    auto Exstar = Xh->element();
    auto Eystar = Xh->element();
    auto Bzstar = Xh->element();
    FSolve( wbdy, dt,Ex, Ey, Bz, Exstar, Eystar, Bzstar );

    Ex.add( dt, Exstar );
    Ey.add( dt, Eystar );
    Bz.add( dt, Bzstar );
    time += dt;
}
template<typename BdyExpr>
void
Diode::EulerModifiedStep( double& time,double dt,
                          BdyExpr& wbdy,
                          element_type& Ex, element_type& Ey,element_type& Bz )
{
    element_type Exstar = Xh->element();
    element_type Eystar = Xh->element();
    element_type Bzstar = Xh->element();
    element_type Exn = Ex;
    element_type Eyn = Ey;
    element_type Bzn = Bz;
    FSolve( wbdy, dt, Ex, Ey, Bz, Exstar, Eystar, Bzstar );

    Exn.add( dt/2.0, Exstar );
    Eyn.add( dt/2.0, Eystar );
    Bzn.add( dt/2.0, Bzstar );
    time += dt/2.0;

    FSolve( wbdy, dt, Exn, Eyn, Bzn, Exstar, Eystar, Bzstar );

    Ex.add( dt, Exstar );
    Ey.add( dt, Eystar );
    Bz.add( dt, Bzstar );
    time += dt/2.0;

}
template<typename BdyExpr>
void
Diode::HeunStep( double& time,double dt,
                 BdyExpr& wbdy,
                 element_type& Ex, element_type& Ey,element_type& Bz )
{
    element_type Exstar = Xh->element();
    element_type Eystar = Xh->element();
    element_type Bzstar = Xh->element();
    element_type Exn = Ex;
    element_type Eyn = Ey;
    element_type Bzn = Bz;
    FSolve( wbdy, dt, Ex, Ey, Bz, Exstar, Eystar, Bzstar );

    Ex.add( dt/2.0, Exstar );
    Ey.add( dt/2.0, Eystar );
    Bz.add( dt/2.0, Bzstar );
    Exn.add( dt, Exstar );
    Eyn.add( dt, Eystar );
    Bzn.add( dt, Bzstar );
    time += dt;

    FSolve( wbdy, dt, Exn, Eyn, Bzn, Exstar, Eystar, Bzstar );

    Ex.add( dt/2.0, Exstar );
    Ey.add( dt/2.0, Eystar );
    Bz.add( dt/2.0, Bzstar );

}
template<typename BdyExpr>
void
Diode::RK4Step( double& time,double dt,
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
    FSolve( wbdy, dt, Ex, Ey, Bz, Exk, Eyk, Bzk );

    Ex.add( dt/6.0, Exk );
    Ey.add( dt/6.0, Eyk );
    Bz.add( dt/6.0, Bzk );
    Exstar=Exn;
    Eystar=Eyn;
    Bzstar=Bzn;
    Exstar.add( dt/2.0, Exk );
    Eystar.add( dt/2.0, Eyk );
    Bzstar.add( dt/2.0, Bzk );

    //k2
    time += dt/2.0;
    FSolve( wbdy, dt, Exstar, Eystar, Bzstar, Exk, Eyk, Bzk );

    Ex.add( dt/3.0, Exk );
    Ey.add( dt/3.0, Eyk );
    Bz.add( dt/3.0, Bzk );
    Exstar=Exn;
    Eystar=Eyn;
    Bzstar=Bzn;
    Exstar.add( dt/2.0, Exk );
    Eystar.add( dt/2.0, Eyk );
    Bzstar.add( dt/2.0, Bzk );

    //k3
    //FSolve(wbdy, Exstar, Eystar, Bzstar, Exk, Eyk, Bzk);

    Ex.add( dt/3.0, Exk );
    Ey.add( dt/3.0, Eyk );
    Bz.add( dt/3.0, Bzk );
    Exstar=Exn;
    Eystar=Eyn;
    Bzstar=Bzn;
    Exstar.add( dt, Exk );
    Eystar.add( dt, Eyk );
    Bzstar.add( dt, Bzk );

    //k4
    time += dt/2.0;
    FSolve( wbdy, dt, Exstar, Eystar, Bzstar, Exk, Eyk, Bzk );

    Ex.add( dt/6.0, Exk );
    Ey.add( dt/6.0, Eyk );
    Bz.add( dt/6.0, Bzk );

}

void
Diode::run( const double* X, unsigned long P, double* Y, unsigned long N )
{
    if ( !this->vm().count( "nochdir" ) )
        Environment::changeRepository( boost::format( "examples/maxwell/%1%/%2%/P%3%G%4%/h_%5%/" )
                                       % this->about().appName()
                                       % convex
                                       % Order
                                       % OrderGeo
                                       % X[0] );

    mesh = createGMSHMesh( _mesh=new mesh_type,
                           _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                           _desc=diodegeo( X[0],OrderGeo,convex ) );

    /**
     * The function space and some associated elements(functions) are then defined
     */
    /** \code */
    int counter=0;
    int exportf = this->vm()["exportfreq"].as<int>();
    Xh = space_type::New( mesh );
    auto Xhc = c_space_type::New( mesh );
    auto Ex = Xh->element();
    auto Ey = Xh->element();
    auto Bz = Xh->element();
    auto Exc = Xhc->element();
    auto Eyc = Xhc->element();
    auto Bzc = Xhc->element();
    auto u = Xh->element();
    auto v = Xh->element();

    using namespace Feel::vf;
    double dt=std::min( 0.1*meshSize/( Order+1 ),timeStep );
    //double dt = timeStep;

    double pi=M_PI;
    double k=2*pi; //=0;
    double theta=this->vm()["theta"].as<double>();
    //theta=0;
    double vu=cos( theta );
    double vv=sin( theta );
    double time = 0;
    auto c=cos( k * ( vu * Px() + vv * Py() - cst_ref( time ) ) + cst( pi/2.0 ) );
    auto Ex_exact = -vv*c;
    auto Ey_exact = vu*c;
    auto Bz_exact = c;
    auto w_exact = vec( Ex_exact, Ey_exact, Bz_exact );
    F_Ex = M_backend->newVector( Xh );
    F_Ey = M_backend->newVector( Xh );
    F_Bz = M_backend->newVector( Xh );

    // left hand side
    D=M_backend->newMatrix( Xh, Xh );
    auto a = form2( _test=Xh, _trial=Xh, _matrix=D, _init=true );
    a = integrate( elements( mesh ), idt( Ex )*id( u ) );
    //D->printMatlab("mass.m");
    auto backend = backend_type::build( soption("backend") );
    auto exporter( export_type::New( this->vm(),
                                     ( boost::format( "%1%-%2%" )
                                       % this->about().appName()
                                       % convex ).str() ) );
    auto L2ProjDisc = projector( Xh, Xh );

    if ( this->vm()["initfields"].as<int>() == 1 )
    {
        Ex = L2ProjDisc->project( Ex_exact );
        Ey = L2ProjDisc->project( Ey_exact );
        Bz = L2ProjDisc->project( Bz_exact );
    }

    auto L2Proj = projector( Xhc, Xhc );
    Exc = L2Proj->project( idv( Ex ) );
    Eyc = L2Proj->project( idv( Ey ) );
    Bzc = L2Proj->project( idv( Bz ) );

    exporter->step( time )->setMesh( mesh );

    exporter->step( time )->add( "Exc", L2Proj->project( idv( Ex ) ) );
    exporter->step( time )->add( "Eyc", L2Proj->project( idv( Ey ) ) );
    exporter->step( time )->add( "Bzc", L2Proj->project( idv( Bz ) ) );

    auto w = vec( idv( Ex ),idv( Ey ),idv( Bz ) );
    auto wR = vec( rightfacev( idv( Ex ) ),rightfacev( idv( Ey ) ),rightfacev( idv( Bz ) ) );
    auto wL = vec( leftfacev( idv( Ex ) ),leftfacev( idv( Ey ) ),leftfacev( idv( Bz ) ) );

    std::cout<<(int)rkmethod<<std::endl;

    switch ( rkmethod )
    {
    case EULER_EXPLICIT:
        std::cout << "Euler explicit" << std::endl;
        break;

    case EULER_MODIFIED:
        std::cout << "Euler modified" << std::endl;
        break;

    case HEUN:
        std::cout << "Heun" << std::endl;
        break;

    case RK4:
        std::cout << "RK4" << std::endl;
        break;
    }

    exporter->save();
    std::cout << "Saved initial/exact solution\n";

    for ( time = 0; time <= Tfinal; )
    {
        std::cout << "============================================================" << std::endl;
        std::cout << "time = " << time << "s, dt=" << dt << ", final time=" << Tfinal << ",hsize = "<< meshSize <<", method : " << (int)rkmethod <<std::endl;

        switch ( rkmethod )
        {
        case EULER_EXPLICIT:
            std::cout<<"euler"<<std::endl;
            EulerStep( time, dt, w_exact, Ex, Ey, Bz );
            break;

        case EULER_MODIFIED:
            std::cout<<"euler modif"<<std::endl;
            EulerModifiedStep( time, dt, w_exact, Ex, Ey, Bz );
            break;

        case HEUN:
            std::cout<<"heun"<<std::endl;
            HeunStep( time, dt, w_exact, Ex, Ey, Bz );
            break;

        case RK4:
            std::cout<<"rk4"<<std::endl;
            RK4Step( time, dt, w_exact, Ex, Ey, Bz );
            break;
        }

        /*
        std::cout << "||exact-Ex||_2 = "
                  << integrate(elements(mesh),
                  (idv(Ex)-Ex_exact)*(idv(Ex)-Ex_exact) ).evaluate().norm()
                  << std::endl;
        std::cout << "||exact-Ey||_2 = "
                  << integrate(elements(mesh),
                  (idv(Ey)-Ey_exact)*(idv(Ey)-Ey_exact) ).evaluate().norm()
                  << std::endl;
        std::cout << "||exact-Bz||_2 = "
                  << integrate(elements(mesh),
                  (idv(Bz)-Bz_exact)*(idv(Bz)-Bz_exact) ).evaluate().norm()
                  << std::endl;
        */
        // save
        counter++;

        if ( counter % exportf == 0 )
        {
            std::cout<<"exporting results"<<std::endl;
            exporter->step( time )->setMesh( mesh );

            Exc = L2Proj->project( idv( Ex ) );
            Eyc = L2Proj->project( idv( Ey ) );
            Bzc = L2Proj->project( idv( Bz ) );
            exporter->step( time )->add( "Exc", Exc );
            exporter->step( time )->add( "Eyc", Eyc );
            exporter->step( time )->add( "Bzc", Bzc );

            exporter->save();
        }
    }



} // Diode::run

/**
 * main function: entry point of the program
 */
int
main( int argc, char** argv )
{
    using namespace Feel;

    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=makeAbout() );

    Application app;
    app.add( new Diode );
    app.run();
}






