/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2009-01-04
  Copyright (C) 2009 Christophe Prud'homme
  Copyright (C) 2009-2010 Université Joseph Fourier (Grenoble I)

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
   \file stokes_Poiseuille.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2009-01-04
 */
#include <feel/options.hpp>
#include <feel/feelcore/application.hpp>

#include <feel/feelalg/backend.hpp>

#include <feel/feeldiscr/functionspace.hpp>

#include <feel/feelpoly/im.hpp>

#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelpoly/lagrange.hpp>
#include <feel/feelpoly/crouzeixraviart.hpp>



#include <feel/feelmesh/elements.hpp>

#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/geotool.hpp>

/**
 * This routine returns the list of options using the
 * boost::program_options library. The data returned is typically used
 * as an argument of a Feel::Application subclass.
 *
 * \return the list of options
 */
inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description stokesDirichletDirichletOptions( "Stokes_Poiseuille options" );
    stokesDirichletDirichletOptions.add_options()
    ( "penal", Feel::po::value<double>()->default_value( 0.5 ), "penalisation parameter" )
    ( "f", Feel::po::value<double>()->default_value( 0 ), "forcing term" )
    ( "mu", Feel::po::value<double>()->default_value( 1.0 ), "reaction coefficient component" )
    ( "hsize", Feel::po::value<double>()->default_value( 0.1 ), "first h value to start convergence" )
    ( "bctype", Feel::po::value<int>()->default_value( 0 ), "0 = strong Dirichlet, 1 = weak Dirichlet" )
    ( "bccoeff", Feel::po::value<double>()->default_value( 100.0 ), "coeff for weak Dirichlet conditions" )
    ( "export-matlab", "export matrix and vectors in matlab" )
    ;
    return stokesDirichletDirichletOptions.add( Feel::feel_options() ) ;
}


/**
 * This routine defines some information about the application like
 * authors, version, or name of the application. The data returned is
 * typically used as an argument of a Feel::Application subclass.
 *
 * \return some data about the application.
 */
inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "stokes_Poiseuille" ,
                           "stokes_Poiseuille" ,
                           "0.1",
                           "Stokes equation on simplices or simplex products",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2009-2012 Université de Grenoble 1 (Joseph Fourier)" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "" );
    return about;

}


namespace Feel
{
using namespace Feel::vf;

/**
 * \class Stokes class
 * \brief solves the stokes equations
 *
 */
template<int POrder=1,int GeoOrder=1>
class Stokes_Poiseuille
    :
public Application
{
    typedef Application super;
public:


    typedef double value_type;

    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    /*matrix*/
    typedef backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef backend_type::vector_ptrtype vector_ptrtype;

    /*mesh*/
    /* *********************************************************************
     * STOKESMESHTYPE : mesh type
     * 1 : hypercube,2D, default
     * 2 : Cylinder, 3D
     *************************************************************************/
#ifndef STOKESPRESSMESHTYPE
#define  STOKESPRESSMESHTYPE 1
#endif
#if (STOKESPRESSMESHTYPE == 1)
    typedef Simplex<2> convex_type;
#elif (STOKESPRESSMESHTYPE == 2)

    typedef Simplex<3,GeoOrder> convex_type;
#endif

#ifndef FEELPP_COND
#define  FEELPP_COND DD
#endif

    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    /*basis*/

    typedef Lagrange<POrder+1, Vectorial> basis_u_type;
    typedef Lagrange<POrder, Scalar> basis_p_type;
    typedef Lagrange<0, Scalar> basis_l_type;

    // use lagrange multipliers to ensure zero mean pressure
    #if defined( FEELPP_USE_LM )
    typedef bases<basis_u_type,basis_p_type, basis_l_type> basis_type;
    #else
    typedef bases<basis_u_type,basis_p_type> basis_type;
    #endif

    typedef bases<basis_u_type> basis_type_U;
    typedef FunctionSpace<mesh_type, basis_type_U> space_type_U;
    typedef boost::shared_ptr<space_type_U> space_ptrtype_U;

    /*space*/
    //# marker2 #
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef boost::shared_ptr<space_type> space_ptrtype;
    //# endmarker2 #

    /* functions */
    //# marker3 #
    typedef typename space_type::element_type element_type;
    //# endmarker3 #

    /* export */
    typedef Exporter<mesh_type> export_type;

    FEELPP_DONT_INLINE
    Stokes_Poiseuille( );

    // init mesh and space
    FEELPP_DONT_INLINE
    void init();

    /**
     * run the convergence test
     */
    FEELPP_DONT_INLINE
    void run();

private:


    /**
     * export results to ensight format (enabled by  --export cmd line options)
     */
    template<typename ExprUExact, typename ExprPExact>
    FEELPP_DONT_INLINE
    void exportResults( ExprUExact uexact, ExprPExact pexact,
                        element_type& u, element_type& v );

private:

    backend_ptrtype M_backend;
    double meshSize;

    double mu;
    double penalbc;

    mesh_ptrtype mesh;
    space_ptrtype Xh;
    space_ptrtype_U Uh;
    sparse_matrix_ptrtype M,D;
    vector_ptrtype F;

    boost::shared_ptr<export_type> exporter;
}; // Stokes

template<int POrder, int GeoOrder>
Stokes_Poiseuille<POrder,GeoOrder>::Stokes_Poiseuille( )
    :
    super( ),
    M_backend( backend_type::build( soption("backend") ) ),
    meshSize( doption("hsize") ),
    mu( this->vm()["mu"].template as<value_type>() ),
    penalbc( this->vm()["bccoeff"].template as<value_type>() ),
    exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) )
{

}

template<int POrder, int GeoOrder>
void
Stokes_Poiseuille<POrder,GeoOrder>::init()
{
    if ( this->vm().count( "help" ) )
    {
        std::cout << this->optionsDescription() << "\n";
        return;
    }

    if ( this->vm().count( "nochdir" ) == false )
        this->changeRepository( boost::format( "benchmarks/%1%/%2%/Part%6%/P%3%P%4%/h_%5%/" )
                                % this->about().appName()
                                % convex_type::name()
                                % basis_u_type::nOrder % basis_p_type::nOrder
                                % doption("hsize")
                                % Environment::numberOfProcessors() );

#if (STOKESPRESSMESHTYPE == 1)
    //********************** Rectangle ***************************************
    mesh = createGMSHMesh( _mesh=new mesh_type,
                            _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                            _desc=domain( _name= ( boost::format( "%1%-%2%-%3%" ) % "hypercube" % convex_type().dimension() % 1 ).str() ,
                                          _shape="hypercube",
                                          _dim=convex_type().dimension(),
                                          _h=meshSize ) );

     mesh->addMarkerName( "inlet", 1, 1 );
     mesh->addMarkerName( "outlet", 3, 1 );
     mesh->addMarkerName( "wall1", 2, 1 );
     mesh->addMarkerName( "wall2", 4, 1 );
    std::cout << "number of elements of 2D: " << mesh->numElements() << "\n";
    //************************ Cylindre ****************************************
#elif (STOKESPRESSMESHTYPE == 2 )
    GeoTool::Node Centre(0,0,0);
    GeoTool::Node Rayon( 1);
    GeoTool::Node Dir(1,0,0);
    GeoTool::Node Lg(5,0,0);
    GeoTool::Cylindre C( meshSize,"Cyl",Centre,Dir,Rayon,Lg);
    C.setMarker(_type="surface",_name="inlet",_marker1=true);
    C.setMarker(_type="surface",_name="outlet",_marker2=true);
    C.setMarker(_type="surface",_name="wall",_marker3=true);
    C.setMarker(_type="volume",_name="OmegaFluide",_markerAll=true);

    mesh = C.createMesh(_mesh= new mesh_type,
                        _name="mymesh",
                        _partitions=Environment::worldComm().localSize(),
                        _worldcomm=Environment::worldComm() );
    std::cout << "number of elements : " << mesh->numElements() << "\n";
#endif
    //*************************************************************************





    Xh = space_type::New( mesh );
    Uh=space_type_U::New(mesh);
    F = M_backend->newVector( Xh );
    D =  M_backend->newMatrix( Xh, Xh );
}
template<int POrder, int GeoOrder>
void
Stokes_Poiseuille<POrder,GeoOrder>::run()
{
    this->init();

    auto U = Xh->element( "(u,p)" );
    auto V = Xh->element( "(u,q)" );
    auto u = U.template element<0>( "u" );
    auto v = V.template element<0>( "u" );
    auto p = U.template element<1>( "p" );
    auto q = V.template element<1>( "p" );
    #if defined( FEELPP_USE_LM )
    auto lambda = U.template element<2>();
    auto nu = V.template element<2>();
    #endif
    //# endmarker4 #

    LOG(INFO) << "Data Summary:\n";
    LOG(INFO) << "   hsize = " << meshSize << "\n";
    LOG(INFO) << "  export = " << this->vm().count( "export" ) << "\n";
    LOG(INFO) << "      mu = " << mu << "\n";
    LOG(INFO) << " bccoeff = " << penalbc << "\n";




    //# marker5 #
    auto deft = sym(gradt( u ));
    auto def = sym(grad( v ));
    //# endmarker5 #

    //# marker6 #
    // total stress tensor (trial)
    auto SigmaNt = -idt( p )*N()+mu*deft*N();

    // total stress tensor (test)
    auto SigmaN = -id( p )*N()+mu*def*N();
    //# endmarker6 #


    auto r=1;
    auto L=5;
    auto P_inlet=1;
    auto P_outlet=0.;
#if (STOKESPRESSMESHTYPE == 1)
    //*********************  solution exacte (profile de poiseuille)2D *****************
    auto u_exact=vec((1-(Py()*Py())/(r*r))*(P_inlet-P_outlet)/(2*mu*L),cst(0.) );

    auto p_exact=(P_outlet-P_inlet)*Px()/L + P_inlet;

    auto f=vec(cst(0.) , cst(0.) );
    //*************************************************************************************
#elif  (STOKESPRESSMESHTYPE == 2)
    //*********************  solution exacte (profile de poiseuille)3D **********************
    auto u_exact=vec(  (P_inlet-P_outlet)*r*r*(1-(Py()*Py()+Pz()*Pz())/(r*r))/(4*L) , cst(0.) , cst(0.) );

    auto p_exact=(-Px()*(P_inlet-P_outlet)/L + P_inlet);

    auto f=vec(cst(0.) , cst(0.), cst(0.) );
    //*************************************************************************************
#endif


    double taille = 4*math::atan(1.)*r*r*L;
    double mean_p_exact = integrate( elements( mesh ),  p_exact ).evaluate()(0,0) /taille;
    std::cout << "[stokes] mean(p_exact)=" << mean_p_exact << "\n";

    //***************************CL: Dirichlet-Neumann**********************************
 #if(FEELPP_COND==DN)
    auto u_ex_proj=vf::project(Uh,elements(mesh),u_exact);
    auto D_ex = sym(gradv(u_ex_proj ));

    // right hand side
    auto stokes_rhs = form1( _test=Xh, _vector=F );
    stokes_rhs += integrate( elements( mesh ),inner( f,id( v ) ) );
    LOG(INFO) << "[stokes] vector local assembly done\n";
    stokes_rhs += integrate( markedfaces( mesh,"outlet" ), inner( (-p_exact)*N(),id(v) ) );
    stokes_rhs += integrate( markedfaces( mesh,"outlet" ), inner( (2*D_ex)*N(),id(v) ) );

    // left hand side
    auto stokes = form2( _test=Xh, _trial=Xh, _matrix=D );
    mpi::timer chrono;
    stokes += integrate( elements( mesh ),2*mu*inner( deft,def ) );
    std::cout << "mu*inner(deft,def): " << chrono.elapsed() << "\n";
    chrono.restart();
    stokes +=integrate( elements( mesh ), - div( v )*idt( p ) + divt( u )*id( q ));
    std::cout << "(u,p): " << chrono.elapsed() << "\n";
    chrono.restart();

    stokes+=on( _range=markedfaces(mesh, "inlet"), _element=u,_rhs=stokes_rhs,
                _expr=u_exact );
    #if (STOKESPRESSMESHTYPE == 2)
    stokes+=on( _range=markedfaces(mesh, "wall"), _element=u,_rhs=stokes_rhs,
                _expr=u_exact );
    #elif (STOKESPRESSMESHTYPE == 1)
    stokes+=on( _range=markedfaces(mesh, "wall1"), _element=u,_rhs=stokes_rhs,
                 _expr=u_exact );
    stokes+=on( _range=markedfaces(mesh, "wall2"), _element=u,_rhs=stokes_rhs,
                 _expr=u_exact );
    #endif

    std::cout << "bc: " << chrono.elapsed() << "\n";
    chrono.restart();
    stokes.solve( _rhs=stokes_rhs, _solution=U);

//***************************CL: Dirichlet-Dirichlet **********************************
#elif(FEELPP_COND==DD)
    // right hand side
    auto stokes_rhs = form1( _test=Xh, _vector=F );
    stokes_rhs += integrate( elements( mesh ),inner( f,id( v ) ));
    #if defined( FEELPP_USE_LM )
    stokes_rhs += integrate( elements( mesh ),id(nu)* p_exact );   //add the mean of the exact pressure
    #endif
    //************
    LOG(INFO) << "[stokes] vector local assembly done\n";

    // Construction of the left hand side
    //# marker7 #
    auto stokes = form2( _test=Xh, _trial=Xh, _matrix=D );
    boost::timer chrono;
    stokes += integrate( elements( mesh ), 2*mu*inner( deft,def ) );
    std::cout << "mu*inner(deft,def): " << chrono.elapsed() << "\n";
    chrono.restart();
    stokes +=integrate( elements( mesh ), - div( v )*idt( p ) + divt( u )*id( q ) );
    std::cout << "(u,p): " << chrono.elapsed() << "\n";
    chrono.restart();
    #if defined( FEELPP_USE_LM )
    stokes +=integrate( elements( mesh ), id( q )*idt( lambda ) + idt( p )*id( nu ));
    std::cout << "(lambda,p): " << chrono.elapsed() << "\n";
    chrono.restart();
    #endif

    stokes+=on( _range=boundaryfaces(mesh), _element=u,_rhs=stokes_rhs,
                _expr=u_exact );

    std::cout << "bc: " << chrono.elapsed() << "\n";
    chrono.restart();
    stokes.solve( _rhs=stokes_rhs, _solution=U);

//***************************CL: Neumann-Neumann**********************************
#elif (FEELPP_COND==NN)
    auto u_ex_proj=vf::project(Uh,elements(mesh),u_exact);
    auto D_ex = sym(gradv(u_ex_proj ));

    // right hand side
    auto stokes_rhs = form1( _test=Xh, _vector=F );
    stokes_rhs += integrate( elements( mesh ),inner( f,id( v ) ) );
    LOG(INFO) << "[stokes] vector local assembly done\n";
    stokes_rhs += integrate( markedfaces( mesh,"outlet" ), inner( (-p_exact)*N(),id(v) ) );
    stokes_rhs += integrate( markedfaces( mesh,"outlet" ), inner( (2*D_ex)*N(),id(v) ) );
    stokes_rhs += integrate( markedfaces( mesh,"inlet" ), inner( (-p_exact)*N(),id(v) ) );
    stokes_rhs += integrate( markedfaces( mesh,"inlet" ), inner( (2*D_ex)*N(),id(v) ) );

    // left hand side
    auto stokes = form2( _test=Xh, _trial=Xh, _matrix=D );
    mpi::timer chrono;
    stokes += integrate( elements( mesh ),2*mu*inner( deft,def ) );
    std::cout << "mu*inner(deft,def): " << chrono.elapsed) () << "\n";
    chrono.restart();
    stokes +=integrate( elements( mesh ), - div( v )*idt( p ) + divt( u )*id( q ));
    std::cout << "(u,p): " << chrono.elapsed() << "\n";
    chrono.restart();
    #if (STOKESPRESSMESHTYPE == 2)
    stokes+=on( _range=markedfaces(mesh, "wall"), _element=u,_rhs=stokes_rhs,
                _expr=u_exact );
    #elif (STOKESPRESSMESHTYPE == 1)
    stokes+=on( _range=markedfaces(mesh, "wall1"), _element=u,_rhs=stokes_rhs,
                 _expr=u_exact );
    stokes+=on( _range=markedfaces(mesh, "wall2"), _element=u,_rhs=stokes_rhs,
                 _expr=u_exact );
    #endif
    std::cout << "bc: " << chrono.elapsed() << "\n";
    chrono.restart();,
    stokes.solve( _rhs=stokes_rhs, _solution=U);
#endif

size_type nnz = 0 ;
auto nNz = D->graph()->nNz() ;
for ( auto iter = nNz.begin(); iter!=nNz.end(); ++iter )
    nnz += ( *iter ) ;
size_type gnnz=0;
mpi::all_reduce( this->comm(), nnz, gnnz, [] ( size_type x, size_type y ) {return x + y;} );
LOG(INFO) << "[nnz]       number of local nnz: " << nnz << " global: "  << gnnz << "\n";
LOG(INFO) << "[dof]             number of dof: " << Xh->nDof() << "\n";
LOG(INFO) << "[dof]        number of dof/proc: " << Xh->nLocalDof() << "\n";
LOG(INFO) << "[dof]          number of dof(U): " << Xh->template functionSpace<0>()->nDof()  << "\n";
LOG(INFO) << "[dof]     number of dof/proc(U): " << Xh->template functionSpace<0>()->nLocalDof()  << "\n";
LOG(INFO) << "[dof]          number of dof(P): " << Xh->template functionSpace<1>()->nDof()  << "\n";
LOG(INFO) << "[dof]     number of dof/proc(P): " << Xh->template functionSpace<1>()->nLocalDof()  << "\n";

std::cout << "[nnz]       number of local nnz: " << nnz << " global: "  << gnnz << "\n";
std::cout << "[dof]             number of dof: " << Xh->nDof() << "\n";
std::cout << "[dof]        number of dof/proc: " << Xh->nLocalDof() << "\n";
std::cout << "[dof]          number of dof(U): " << Xh->template functionSpace<0>()->nDof()  << "\n";
std::cout << "[dof]     number of dof/proc(U): " << Xh->template functionSpace<0>()->nLocalDof()  << "\n";
std::cout << "[dof]          number of dof(P): " << Xh->template functionSpace<1>()->nDof()  << "\n";
std::cout << "[dof]     number of dof/proc(P): " << Xh->template functionSpace<1>()->nLocalDof()  << "\n";


 this->exportResults( u_exact, p_exact, U, V );
} // Stokes::run


template<int POrder, int GeoOrder>
template<typename ExprUExact, typename ExprPExact>
void
Stokes_Poiseuille<POrder,GeoOrder>::exportResults( ExprUExact u_exact, ExprPExact p_exact,
                       element_type& U, element_type& V )
{
    auto u = U.template element<0>();
    auto p = U.template element<1>();

    auto v = V.template element<0>();
    auto q = V.template element<1>();
    #if defined( FEELPP_USE_LM )
         auto lambda = U.template element<2>();
         auto nu = V.template element<2>();
         LOG(INFO) << "value of the Lagrange multiplier lambda= " << lambda( 0 ) << "\n";
         std::cout << "value of the Lagrange multiplier lambda= " << lambda( 0 ) << "\n";
    #endif

    auto u_exact_proj=vf::project(Uh,elements(mesh),u_exact);

    double u_errorL2 = integrate( elements( u.mesh() ), trans( idv( u )-u_exact )*( idv( u )-u_exact ) ).evaluate()( 0, 0 );
    std::cout << "||u_error||_2 = " << math::sqrt( u_errorL2 ) << "\n";
    LOG(INFO) <<"||u_error||_2 = " << math::sqrt( u_errorL2 ) << "\n";

    double meas = integrate( elements( u.mesh() ), cst( 1.0 )).evaluate()( 0, 0 );
#if (STOKESPRESSMESHTYPE ==2)
    LOG(INFO) << "[stokes] measure(Omega)=" << meas << " (should be equal to "<< 4*math::atan(1.)*5 << ")\n";
    std::cout << "[stokes] measure(Omega)=" << meas << " (should be equal to "<< 4*math::atan(1.)*5 << ")\n";
#elif (STOKESPRESSMESHTYPE == 1)
    LOG(INFO) << "[stokes] measure(Omega)=" << meas << " (should be equal to  1)\n";
    std::cout << "[stokes] measure(Omega)=" << meas << " (should be equal to  1)\n";
#endif

    double mean_p = integrate( elements( u.mesh() ), idv( p ) ).evaluate()( 0, 0 )/meas;
    LOG(INFO) << "[stokes] mean(p)=" << mean_p << "\n";

    std::cout << "[stokes] mean(p)=" << mean_p << "\n";


    // double p_errorL2 = integrate( elements( u.mesh() ), ( idv( p )+mean_p_exact - p_exact )*( idv( p )+mean_p_exact-p_exact ), _quad=_Q<QuadOrder>() ).evaluate()( 0, 0 );
    double p_errorL2 = integrate( elements( u.mesh() ), ( idv( p ) - p_exact )*( idv( p )-p_exact )).evaluate()( 0, 0 );
    std::cout << "||p_error||_2 = " << math::sqrt( p_errorL2 ) << "\n";
    LOG(INFO) <<"||p_error||_2 = " << math::sqrt( p_errorL2 ) << "\n";
    LOG(INFO) << "[stokes] solve for D done\n";



    double u_errorH1 = integrate( elements( u.mesh() ),  trans( idv( u )-u_exact )*( idv( u )-u_exact )).evaluate()( 0, 0 ) +  integrate( elements( u.mesh() ),  trans( gradv( u ) -  gradv ( u_exact_proj ) )*( gradv( u ) -  gradv ( u_exact_proj )) ).evaluate()( 0, 0 );
    double H1_u=math::sqrt( u_errorH1 );
    std::cout << "||u_errorH1||_2 = " <<H1_u<< "\n";
    LOG(INFO) << "||u_errorH1||_2 = " <<H1_u<< "\n";

    double mean_div_u = integrate( elements( u.mesh() ), divv( u ) ).evaluate()( 0, 0 );
    LOG(INFO) << "[stokes] mean_div(u)=" << mean_div_u << "\n";
    std::cout << "[stokes] mean_div(u)=" << mean_div_u << "\n";

    double div_u_error_L2 = integrate( elements( u.mesh() ), divv( u )*divv( u ) ).evaluate()( 0, 0 );
    LOG(INFO) << "[stokes] ||div(u)||_2=" << math::sqrt( div_u_error_L2 ) << "\n";

    std::cout << "[stokes] ||div(u)||=" << math::sqrt( div_u_error_L2 ) << "\n";

    v = vf::project( u.functionSpace(), elements( u.mesh() ), u_exact );
    q = vf::project( p.functionSpace(), elements( p.mesh() ), p_exact );

#if (STOKESPRESSMESHTYPE ==2)
    auto deff = sym(gradv(u));
    auto SigmaNN =(-idv(p)*vf::N()+2*mu*deff*vf::N());
    auto FappIn = integrate(markedfaces( mesh,"inlet" ) , SigmaNN).evaluate();
    auto Fapp1In = FappIn(0,0); //pour la prémière composante
    auto Fapp2In = FappIn(1,0); //pour la seconde
    auto Fapp3In = FappIn(2,0);
    std::cout << "Fapp1In = "<<math::abs(pi- FappIn(0,0)) << "\n" ;
    std::cout << "Fapp2In = "<<math::abs( -FappIn(1,0)) << "\n" ;
    std::cout << "Fapp3In = "<<math::abs( -FappIn(2,0)) << "\n" ;
    LOG(INFO) << "Fapp1In = "<<math::abs( pi-FappIn(0,0)) << "\n" ;
    LOG(INFO) << "Fapp2In = "<<math::abs( -FappIn(1,0)) << "\n" ;
    LOG(INFO) << "Fapp3In = "<<math::abs( -FappIn(2,0)) << "\n" ;

    auto FappOut = integrate(markedfaces( mesh,"outlet" ) , SigmaNN).evaluate();
    auto Fapp1Out = FappOut(0,0); //pour la prémière composante
    auto Fapp2Out = FappOut(1,0); //pour la seconde
    auto Fapp3Out = FappOut(2,0);
    std::cout << "Fapp1Out = "<<math::abs( -FappOut(0,0)) << "\n" ;
    std::cout << "Fapp2Out = "<<math::abs( -FappOut(1,0)) << "\n" ;
    std::cout << "Fapp3Out = "<<math::abs( -FappOut(2,0)) << "\n" ;
    LOG(INFO) << "Fapp1Out = "<<math::abs( -FappOut(0,0)) << "\n" ;
    LOG(INFO) << "Fapp2Out = "<<math::abs( -FappOut(1,0)) << "\n" ;
    LOG(INFO) << "Fapp3Out = "<<math::abs( -FappOut(2,0)) << "\n" ;

    auto FappWall = integrate(markedfaces( mesh,"wall" ) , SigmaNN).evaluate();
    auto Fapp1Wall = FappWall(0,0); //pour la prémière composante
    auto Fapp2Wall = FappWall(1,0); //pour la seconde
    auto Fapp3Wall = FappWall(2,0);
    std::cout << "Fapp1Wall = "<< math::abs(-pi-FappWall(0,0)) << "\n" ;
    std::cout << "Fapp2Wall = "<<math::abs( -FappWall(1,0) )<< "\n" ;
    std::cout << "Fapp3Wall = "<< math::abs(-FappWall(2,0)) << "\n" ;
    LOG(INFO) << "Fapp1Wall = "<< math::abs(-pi-FappWall(0,0)) << "\n" ;
    LOG(INFO) << "Fapp2Wall = "<< math::abs(-FappWall(1,0)) << "\n" ;
    LOG(INFO) << "Fapp3Wall = "<<math::abs( -FappWall(2,0)) << "\n" ;

#endif


    if ( exporter->doExport() )
    {
        exporter->step( 0 )->setMesh( U.functionSpace()->mesh() );
        exporter->step( 0 )->addRegions();
        auto v = U.functionSpace()->template functionSpace<0> ()->element();
        v = U.template element<0>();
        exporter->step( 0 )->add( "u", U.template element<0>() );
        //exporter->step( 0 )->add( "ux", v.comp( X ) );
        //exporter->step( 0 )->add( "uy", v.comp( Y ) );
        exporter->step( 0 )->add( "u", U.template element<0>() );
        exporter->step( 0 )->add( "p", U.template element<1>() );
        exporter->step( 0 )->add( "u_exact", V.template element<0>() );
        exporter->step( 0 )->add( "p_exact", V.template element<1>() );
        exporter->save();
    }

} // Stokes::export
} // Feel

int
main( int argc, char** argv )
{

    using namespace Feel;

    Environment env( _argc=argc, _argv=argv,
                     _desc=makeOptions(),
                     _about=about(_name="stokes_Poiseuille",
                                  _author="Christophe Prud'homme",
                                  _email="christophe.prudhomme@feelpp.org") );

    Feel::Stokes_Poiseuille<1,1> Stokes_Poiseuille;
    Stokes_Poiseuille.run();
}











