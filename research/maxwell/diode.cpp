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
    static const uint16_type Order = 2;

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
    typedef Simplex<2,Order> convex_type;
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
    typedef Exporter<mesh_type> export_type;
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
        convex( this->vm()["convex"].as<std::string>() )
    {
    }

    void run();

    void run( const double* X, unsigned long P, double* Y, unsigned long N );

private:

    //! linear algebra backend
    backend_ptrtype M_backend;

    //! mesh characteristic size
    double meshSize;
    double timeStep;
    double Tfinal;
    //! convex of the domain
    std::string convex;
}; // Diode
const uint16_type Diode::Order;

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
void
Diode::run( const double* X, unsigned long P, double* Y, unsigned long N )
{
    if ( !this->vm().count( "nochdir" ) )
        Environment::changeRepository( boost::format( "examples/maxwell/%1%/%2%/P%3%/h_%4%/" )
                                       % this->about().appName()
                                       % convex
                                       % Order
                                       % X[0] );
    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _desc=geo(_filename=(boost::format("diode-%1%.geo")%convex).str(),
                                                  _dim=2,
                                                  _order=Order,
                                                  _h=X[0] ) );

    /**
     * The function space and some associated elements(functions) are then defined
     */
    /** \code */
    auto Xh = space_type::New( mesh );
    auto Xhc = c_space_type::New( mesh );
    auto Ex = Xh->element();
    auto Ey = Xh->element();
    auto Bz = Xh->element();
    auto Exe = Xhc->element();
    auto Eye = Xhc->element();
    auto Bze = Xhc->element();
    auto Exc = Xhc->element();
    auto Eyc = Xhc->element();
    auto Bzc = Xhc->element();
    auto u = Xh->element();
    auto v = Xh->element();

    using namespace Feel::vf;
    double dt=std::min(0.1*meshSize/Order,timeStep);

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
    auto Anp_1 = vec( +Ny() * Ny() / 0.2e1, -Nx() * Ny() / 0.2e1, -Ny() / 0.2e1 );
    auto Anp_2 = vec(-Nx() * Ny() / 0.2e1, Nx() * Nx() / 0.2e1, Nx() / 0.2e1 );
    auto Anp_3 = vec( -Ny() / 0.2e1, Nx() / 0.2e1, cst(0.1e1 / 0.2e1) );
    auto Anm_1 = vec(  -Ny() * Ny() / 0.2e1, Nx() * Ny() / 0.2e1, -Ny() / 0.2e1 );
    auto Anm_2 = vec( Nx() * Ny() / 0.2e1, -Nx() * Nx() / 0.2e1, Nx() / 0.2e1);
    auto Anm_3 = vec( -Ny() / 0.2e1, Nx() / 0.2e1, cst(-0.1e1 / 0.2e1) );
    auto F_Ex = M_backend->newVector( Xh );
    auto lEx=form1( _test=Xh, _vector=F_Ex, _init=true );
    auto F_Ey = M_backend->newVector( Xh );
    auto lEy=form1( _test=Xh, _vector=F_Ey, _init=true );
    auto F_Bz = M_backend->newVector( Xh );
    auto lBz=form1( _test=Xh, _vector=F_Bz, _init=true );

    // left hand side
    auto D=M_backend->newMatrix( Xh, Xh );
    auto a = form2( _test=Xh, _trial=Xh, _matrix=D, _init=true );
    a = integrate( elements(mesh), idt(Ex)*id(u)/dt );
    D->printMatlab("mass.m");
    auto backend = backend_type::build( this->vm() );
    auto exporter( export_type::New( this->vm(),
                                     (boost::format( "%1%-%2%" )
                                      % this->about().appName()
                                      % convex).str() ) );
    Ex = vf::project( Xh, elements(mesh), Ex_exact );
    Ey = vf::project( Xh, elements(mesh), Ey_exact );
    Bz = vf::project( Xh, elements(mesh), Bz_exact );

    auto L2Proj = projector( Xh, Xhc );
    Exc = L2Proj->project( idv(Ex) );
    Eyc = L2Proj->project( idv(Ey) );
    Bzc = L2Proj->project( idv(Bz) );

    Exe = vf::project( Xhc, elements(mesh), Ex_exact );
    Eye = vf::project( Xhc, elements(mesh), Ey_exact );
    Bze = vf::project( Xhc, elements(mesh), Bz_exact );

    exporter->step(time)->setMesh( mesh );
    exporter->step(time)->add( "Ex", Ex );
    exporter->step(time)->add( "Ey", Ey );
    exporter->step(time)->add( "Bz", Bz );
    exporter->step(time)->add( "Exc", Exc );
    exporter->step(time)->add( "Eyc", Eyc );
    exporter->step(time)->add( "Bzc", Bzc );
    exporter->step(time)->add( "ExExact", Exe );
    exporter->step(time)->add( "EyExact", Eye );
    exporter->step(time)->add( "BzExact", Bze );
    exporter->save();
    std::cout << "Saved initial/exact solution\n";
    for( time = dt; time <= Tfinal; time += dt )
    {
        std::cout << "============================================================" << std::endl;
        std::cout << "time = " << time << "s, dt=" << dt << ", final time=" << Tfinal << std::endl;
        auto w = vec(idv(Ex),idv(Ey),idv(Bz));
        auto wR = vec(rightfacev(idv(Ex)),rightfacev(idv(Ey)),rightfacev(idv(Bz)));
        auto wL = vec(leftfacev(idv(Ex)),leftfacev(idv(Ey)),leftfacev(idv(Bz)));
        // update right hand side
        lEx = integrate( elements( mesh ), idv(Ex)*id(u)/dt  -dyv(Bz)*id(u));
        lEx += integrate( internalfaces(mesh),
                          trans(Anm_1)*(wR-wL)*leftface(id(u)) +
                          trans(Anp_1)*(wR-wL)*rightface(id(u)) );
        lEx += integrate( boundaryfaces(mesh), trans(Anm_1)*(w_exact-w)*id(u) );
        lEy = integrate( elements( mesh ), idv(Ey)*id(u)/dt + dxv(Bz)*id(u));
        lEy += integrate( internalfaces(mesh),
                          trans(Anm_2)*(wR-wL)*leftface(id(u)) +
                          trans(Anp_2)*(wR-wL)*rightface(id(u)) );
        lEy += integrate( boundaryfaces(mesh), trans(Anm_2)*(w_exact-w)*id(u) );
        lBz = integrate( elements( mesh ), idv(Bz)*id(u)/dt+(dxv(Ey)-dyv(Ex))*id(u));
        lBz += integrate( internalfaces(mesh),
                          trans(Anm_3)*(wR-wL)*leftface(id(u)) +
                          trans(Anp_3)*(wR-wL)*rightface(id(u)) );
        lBz += integrate( boundaryfaces(mesh), trans(Anm_3)*(w_exact-w)*id(u) );

        backend->solve( _matrix=D, _solution=Ex, _rhs=F_Ex  );
        backend->solve( _matrix=D, _solution=Ey, _rhs=F_Ey  );
        backend->solve( _matrix=D, _solution=Bz, _rhs=F_Bz  );

        std::cout << "||exact-Ex||_2" << integrate(elements(mesh), (idv(Ex)-Ex_exact)*(idv(Ex)-Ex_exact) ).evaluate().norm() << std::endl;
        std::cout << "||exact-Ey||_2" << integrate(elements(mesh), (idv(Ey)-Ey_exact)*(idv(Ey)-Ey_exact) ).evaluate().norm() << std::endl;
        std::cout << "||exact-Bz||_2" << integrate(elements(mesh), (idv(Bz)-Bz_exact)*(idv(Bz)-Bz_exact) ).evaluate().norm() << std::endl;

        // save
        exporter->step(time)->setMesh( mesh );
        exporter->step(time)->add( "Ex", Ex );
        exporter->step(time)->add( "Ey", Ey );
        exporter->step(time)->add( "Bz", Bz );

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






