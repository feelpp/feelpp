/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-01-04

  Copyright (C) 2008 Christophe Prud'homme
  Copyright (C) 2008-2010 Université Joseph Fourier (Grenoble I)

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
   \file stokes.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-01-04
 */
#include <feel/options.hpp>
#include <feel/feelcore/application.hpp>
#include <feel/feelcore/application.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelpoly/im.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/gmshhypercubedomain.hpp>
#include <feel/feelpoly/polynomialset.hpp>

#include <feel/feelalg/backend.hpp>

#include <feel/feelmesh/elements.hpp>

#include <feel/feelvf/vf.hpp>

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description stokesoptions( "Stokes options" );
    stokesoptions.add_options()
    ( "penal", Feel::po::value<double>()->default_value( 0.5 ), "penalisation parameter" )
    ( "f", Feel::po::value<double>()->default_value( 0 ), "forcing term" )
    ( "stab", Feel::po::value<bool>()->default_value( true ), "true to enable stabilisation, false otherwise" )
    ( "bx", Feel::po::value<double>()->default_value( 1.0 ), "convection X component" )
    ( "by", Feel::po::value<double>()->default_value( 0.0 ), "convection Y component" )
    ( "bz", Feel::po::value<double>()->default_value( 0.0 ), "convection Z component" )
    ( "mu", Feel::po::value<double>()->default_value( 1.0 ), "reaction coefficient component" )
    ( "epsilon", Feel::po::value<double>()->default_value( 1.0 ), "diffusion coefficient" )
    ( "hsize", Feel::po::value<double>()->default_value( 0.5 ), "first h value to start convergence" )
    ( "bctype", Feel::po::value<int>()->default_value( 0 ), "0 = strong Dirichlet, 1 = weak Dirichlet" )
    ( "bccoeff", Feel::po::value<double>()->default_value( 100.0 ), "coeff for weak Dirichlet conditions" )
    ( "export", "export results(ensight, data file(1D)" )
    ( "export-matlab", "export matrix and vectors in matlab" )
    ;
    return stokesoptions.add( Feel::feel_options() ) ;
}
inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "stokes" ,
                           "stokes" ,
                           "0.2",
                           "Stokes equation on simplices or simplex products",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2007,2010 University Joseph Fourier Grenoble 1" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;

}


namespace Feel
{
/**
 * Stokes in a cavity
 *
 */
class Stokes
    :
public Application
{
    typedef Application super;
public:

    // -- TYPEDEFS --
    static const int Dim = 2;
    static const int Order = 2;
    typedef double value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*mesh*/
    typedef Simplex<2> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    /*basis*/
    typedef Lagrange<Order, Vectorial> basis_u_type;
    typedef Lagrange<Order-1, Scalar> basis_p_type;
    typedef Lagrange<0, Scalar> basis_l_type;
    typedef bases<basis_u_type, basis_p_type, basis_l_type> basis_type;

    /*space*/
    typedef FunctionSpace<mesh_type, basis_type, value_type> space_type;

    /* export */
    typedef Exporter<mesh_type> export_type;

    Stokes( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        M_backend( backend_type::build( this->vm() ) ),
        meshSize( this->vm()["hsize"].as<double>() ),
        bcCoeff( this->vm()["bccoeff"].as<double>() ),
        exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) ),
        timers(),
        stats()
    {
        LOG(INFO) << "[Stokes] hsize = " << meshSize << "\n";
        LOG(INFO) << "[Stokes] bccoeff = " << bcCoeff << "\n";
        LOG(INFO) << "[Stokes] export = " << this->vm().count( "export" ) << "\n";

        mu = this->vm()["mu"].as<value_type>();
        penalbc = this->vm()["bccoeff"].as<value_type>();
        epsilon = this->vm()["epsilon"].as<value_type>();
        stab = this->vm()["stab"].as<bool>();

    }


    /**
     * create the mesh using mesh size \c meshSize
     */
    mesh_ptrtype createMesh( double meshSize );


    /**
     * run the convergence test
     */
    void run();

private:

    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    void exportResults( space_type::element_type& u );

private:

    backend_ptrtype M_backend;
    double meshSize;
    double bcCoeff;

    double mu;
    double epsilon;
    bool stab;
    double penalbc;

    boost::shared_ptr<export_type> exporter;

    std::map<std::string,std::pair<boost::timer,double> > timers;
    std::map<std::string,double> stats;
}; // Stokes
const int Stokes::Order;
const int Stokes::Dim;

Stokes::mesh_ptrtype
Stokes::createMesh( double meshSize )
{
    timers["mesh"].first.restart();
    mesh_ptrtype mesh( new mesh_type );


    GmshHypercubeDomain td( convex_type::nDim,convex_type::nOrder,convex_type::nRealDim,convex_type::is_hypercube );
    td.setCharacteristicLength( meshSize );
    std::string fname = td.generate( convex_type::name().c_str() );

    ImporterGmsh<mesh_type> import( fname );
    mesh->accept( import );
    timers["mesh"].second = timers["mesh"].first.elapsed();
    LOG(INFO) << "[timer] createMesh(): " << timers["mesh"].second << "\n";
    return mesh;
} // Stokes::createMesh


void
Stokes::run()
{
    if ( this->vm().count( "help" ) )
    {
        std::cout << this->optionsDescription() << "\n";
        return;
    }

    //    int maxIter = 10.0/meshSize;
    using namespace Feel::vf;

    this->changeRepository( boost::format( "%1%/%2%/P%3%/h_%4%/" )
                            % this->about().appName()
                            % convex_type::name()
                            % Order
                            % this->vm()["hsize"].as<double>()
                          );

    /*
     * First we create the mesh
     */
    mesh_ptrtype mesh = createMesh( meshSize );
    stats["nelt"] = mesh->elements().size();

    /*
     * The function space and some associate elements are then defined
     */
    timers["init"].first.restart();
    auto Xh = space_type::New( mesh );
    //Xh->dof()->showMe();
    auto U = Xh->element();
    auto V = Xh->element();
    auto u = U.element<0>();
    auto v = V.element<0>();
    auto p = U.element<1>();
    auto q = V.element<1>();
    auto lambda = U.element<2>();
    auto nu = V.element<2>();

    timers["init"].second = timers["init"].first.elapsed();
    stats["ndof"] = Xh->nDof();

    LOG(INFO) << "Data Summary:\n";
    LOG(INFO) << "      mu = " << mu << "\n";
    LOG(INFO) << " bccoeff = " << penalbc << "\n";
    LOG(INFO) << " epsilon = " << epsilon << "\n";
    LOG(INFO) << "    stab = " << stab << "\n";


    auto F= M_backend->newVector( Xh );
    timers["assembly"].first.restart();
    auto deft = 0.5*( gradt( u )+trans( gradt( u ) ) );
    auto def = 0.5*( grad( v )+trans( grad( v ) ) );
    auto Id = ( mat<Dim,Dim>( cst( 1 ), cst( 0 ), cst( 0 ), cst( 1. ) ) );
    auto SigmaNt = ( -idt( p )*N()+2*mu*deft*N() );
    auto SigmaN = ( -id( p )*N()+2*mu*def*N() );
    auto g= oneX();
    form1( Xh, F, _init=true )  =
        //integrate( elements(mesh), im, trans(vec(cst(0.),cst(0.)))*id(v) ) +
        integrate( markedfaces( mesh,4 ),
                   trans( g )*( -SigmaN+penalbc*id( v )/hFace() ) );

    LOG(INFO) << "[stokes] vector local assembly done\n";
    timers["assembly"].second = timers["assembly"].first.elapsed();
    timers["assembly_F"].second = timers["assembly"].first.elapsed();

    /*
     * Construction of the left hand side
     */
    auto S = M_backend->newMatrix( Xh, Xh );
    form2( Xh, Xh, S, _init=true );

    if ( this->vm().count( "export-matlab" ) )
    {
        S->printMatlab( "S.m" );
    }

    auto D= M_backend->newMatrix( Xh, Xh );
    timers["assembly"].first.restart();

    form2( Xh, Xh, D, _init=true ) = integrate( elements( mesh ), mu*trace( deft*trans( def ) ) );

    form2( Xh, Xh, D ) += integrate( elements( mesh ), - div( v )*idt( p ) + divt( u )*id( q ) );
    form2( Xh, Xh, D ) += integrate( elements( mesh ), id( q )*idt( lambda ) + idt( p )*id( nu ) );
    form2( Xh, Xh, D ) += integrate( boundaryfaces( mesh ), -trans( SigmaNt )*id( v ) );
    form2( Xh, Xh, D ) += integrate( boundaryfaces( mesh ), -trans( SigmaN )*idt( u ) );
    form2( Xh, Xh, D ) += integrate( boundaryfaces( mesh ), +penalbc*trans( idt( u ) )*id( v )/hFace() );

    LOG(INFO) << "[stokes] matrix local assembly done\n";
    LOG(INFO) << "[stokes] vector/matrix global assembly done\n";

    if ( this->vm().count( "export-matlab" ) )
    {
        F->printMatlab( "F.m" );
        D->printMatlab( "D.m" );
    }

#if 0

    if ( M_bctype == 0 )
        form2( Xh, Xh, D ) +=
#if 0
            on( markedfaces( mesh,2 ), u.comp( X ), F, constant( 1. )  )+
            on( markedfaces( mesh,2 ), u.comp( Y ), F, constant( 0. )  )+
#else
            on( markedfaces( mesh,2 ), u, F, oneX()  )+
#endif
            on( markedfaces( mesh,1 ), u, F, cst( 0. )*oneX()+cst( 0. )*oneY() )+
            on( markedfaces( mesh,3 ), u, F, cst( 0. )*oneX()+cst( 0. )*oneY() )+
            on( markedfaces( mesh,4 ), u, F, cst( 0. )*oneX()+cst( 0. )*oneY() );

#endif

    LOG(INFO) << "[stokes] dirichlet condition applied\n";
    timers["assembly"].second += timers["assembly"].first.elapsed();
    timers["assembly_D"].second += timers["assembly"].first.elapsed();

    if ( this->vm().count( "export-matlab" ) )
    {
        F->printMatlab( "F_dir.m" );
        D->printMatlab( "D_dir.m" );
    }

    LOG(INFO) << "[stokes] starting solve for D\n";

    timers["solver"].first.restart();
    backend_type::build( this->vm() )->solve( _matrix=D, _solution=U, _rhs=F );
    timers["solver"].second = timers["solver"].first.elapsed();

    if ( this->vm().count( "export-matlab" ) )
    {

        U.printMatlab( "U.m" );
    }

    LOG(INFO) << "[stokes] solve for D done\n";
    double meas = integrate( elements( mesh ), constant( 1.0 ) ).evaluate()( 0, 0 );
    double mean_p = integrate( elements( mesh ), idv( p ) ).evaluate()( 0, 0 )/meas;
    LOG(INFO) << "[stokes] mean(p)=" << mean_p << "\n";

    this->exportResults( U );

    LOG(INFO) << "[dof]         number of dof: " << Xh->nDof() << "\n";
    LOG(INFO) << "[dof]    number of dof/proc: " << Xh->nLocalDof() << "\n";
    LOG(INFO) << "[dof]      number of dof(U): " << Xh->functionSpace<0>()->nDof()  << "\n";
    LOG(INFO) << "[dof] number of dof/proc(U): " << Xh->functionSpace<0>()->nLocalDof()  << "\n";
    LOG(INFO) << "[dof]      number of dof(P): " << Xh->functionSpace<1>()->nDof()  << "\n";
    LOG(INFO) << "[dof] number of dof/proc(P): " << Xh->functionSpace<1>()->nLocalDof()  << "\n";
    LOG(INFO) << "[timer] run():         init: " << timers["init"].second << "\n";
    LOG(INFO) << "[timer] run():     assembly: " << timers["assembly"].second << "\n";
    LOG(INFO) << "[timer] run():         o D : " << timers["assembly_D"].second << "\n";
    LOG(INFO) << "[timer] run():         o F : " << timers["assembly_F"].second << "\n";
    LOG(INFO) << "[timer] run():         o M : " << timers["assembly_M"].second << "\n";
    LOG(INFO) << "[timer] run():         o L : " << timers["assembly_L"].second << "\n";
    LOG(INFO) << "[timer] run():         o i : " << timers["assembly_evaluate"].second << "\n";
    LOG(INFO) << "[timer] run():       solver: " << timers["solver"].second << "\n";
    LOG(INFO) << "[timer] run():       solver: " << timers["export"].second << "\n";

} // Stokes::run


void
Stokes::exportResults( space_type::element_type& U )
{
    timers["export"].first.restart();

    exporter->step( 1. )->setMesh( U.functionSpace()->mesh() );
    exporter->step( 1. )->add( "u", U.element<0>() );
    exporter->step( 1. )->add( "p", U.element<1>() );
    exporter->save();

    timers["export"].second = timers["export"].first.elapsed();
    LOG(INFO) << "[timer] exportResults(): " << timers["export"].second << "\n";
} // Stokes::export
} // Feel




int
main( int argc, char** argv )
{
    Feel::Environment env( argc, argv );
    /* assertions handling */
    Feel::Assert::setLog( "stokes.assert" );

    /* define and run application */
    Feel::Stokes stokes( argc, argv, makeAbout(), makeOptions() );
    stokes.run();
}





