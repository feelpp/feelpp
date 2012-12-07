/*-*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

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
   \file Stokes_Neumann_Neumann.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2009-01-04
 */
#include <feel/options.hpp>
#include <feel/feelcore/application.hpp>

#include <feel/feelalg/backend.hpp>

#include <feel/feeldiscr/functionspace.hpp>

#include <feel/feelpoly/im.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/gmshhypercubedomain.hpp>
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
    Feel::po::options_description stokes_Neumann_Neumann( "Stokes_Neumann_Neumann options" );
    stokes_Neumann_Neumann.add_options()
    ( "penal", Feel::po::value<double>()->default_value( 0.5 ), "penalisation parameter" )
    ( "f", Feel::po::value<double>()->default_value( 0 ), "forcing term" )
    ( "mu", Feel::po::value<double>()->default_value( 1.0 ), "reaction coefficient component" )
    ( "hsize", Feel::po::value<double>()->default_value( 0.1 ), "first h value to start convergence" )
    ( "bctype", Feel::po::value<int>()->default_value( 0 ), "0 = strong Dirichlet, 1 = weak Dirichlet" )
    ( "bccoeff", Feel::po::value<double>()->default_value( 100.0 ), "coeff for weak Dirichlet conditions" )
    ( "export-matlab", "export matrix and vectors in matlab" )
    ;
    return stokes_Neumann_Neumann.add( Feel::feel_options() ) ;
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
    Feel::AboutData about( "Stokes_Neumann_Neumann" ,
                           "Stokes_Neumann_Neumann" ,
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
class Stokes_Neumann_Neumann
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
     * 3 : read aorta mesh, 3D,!!!!! check path to mesh file !!!!!
     *************************************************************************/
#ifndef STOKESPRESSMESHTYPE
#define  STOKESPRESSMESHTYPE 1
#endif
#if (STOKESPRESSMESHTYPE == 1)
    typedef Simplex<2> convex_type;
#elif (STOKESPRESSMESHTYPE == 2)
    typedef Simplex<3> convex_type;
#elif (STOKESPRESSMESHTYPE == 3)
    typedef Simplex<3> convex_type;
#endif

    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    /*basis*/
    //# marker1 #,
    typedef Lagrange<4, Vectorial> basis_u_type;
    typedef Lagrange<3, Scalar> basis_p_type;
    typedef Lagrange<0, Scalar> basis_l_type;

    // use lagrange multipliers to ensure zero mean pressure
#if defined( FEELPP_USE_LM )
    typedef bases<basis_u_type,basis_p_type, basis_l_type> basis_type;
#else
    typedef bases<basis_u_type,basis_p_type> basis_type;
#endif
    //# endmarker1 #

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
    typedef space_type::element_type element_type;
    //# endmarker3 #

    /* export */
    typedef Exporter<mesh_type> export_type;

    FEELPP_DONT_INLINE
    Stokes_Neumann_Neumann( int argc, char** argv, AboutData const& ad, po::options_description const& od );

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

#if (STOKESPRESSMESHTYPE == 1) ||( STOKESPRESSMESHTYPE == 2 )
    template<typename ExprUExact, typename ExprPExact>
    FEELPP_DONT_INLINE
    void exportResults( ExprUExact uexact, ExprPExact pexact,
                        element_type& u, element_type& v );
#elif (STOKESPRESSMESHTYPE == 3 )
    /* no exact solution know for the aorta case */
    FEELPP_DONT_INLINE
    void exportResults(element_type& u, element_type& v );
#endif

private:

    backend_ptrtype M_backend;
    double meshSize;

    double mu;
    double penalbc;

    mesh_ptrtype mesh;
    space_ptrtype Xh;
    space_ptrtype_U P7;
    sparse_matrix_ptrtype M,D;
    vector_ptrtype F;

    boost::shared_ptr<export_type> exporter;
}; // Stokes


Stokes_Neumann_Neumann::Stokes_Neumann_Neumann( int argc, char** argv, AboutData const& ad, po::options_description const& od )
    :
    super( argc, argv, ad, od ),
    M_backend( backend_type::build( this->vm() ) ),
    meshSize( this->vm()["hsize"].as<double>() ),
    mu( this->vm()["mu"].as<value_type>() ),
    penalbc( this->vm()["bccoeff"].as<value_type>() ),
    exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) )
{

}

void
Stokes_Neumann_Neumann::init()
{
    if ( this->vm().count( "help" ) )
    {
        std::cout << this->optionsDescription() << "\n";
        return;
    }



    if ( this->vm().count( "nochdir" ) == false )
        this->changeRepository( boost::format( "doc/tutorial/%1%/aorta_M2/Part%6%/P%3%P%4%/h_%5%/" )
                                % this->about().appName()
                                % convex_type::name()
                                % basis_u_type::nOrder % basis_p_type::nOrder
                                % this->vm()["hsize"].as<double>()
                                % Environment::numberOfProcessors() );

    auto r=1;
    auto L=5;
    /************************************************************************
     ********************* MESH GENERATION **********************************
     ***********************************************************************/
#if (STOKESPRESSMESHTYPE == 1)
    //********************** Rectangle ***************************************
    mesh = createGMSHMesh( _mesh=new mesh_type,
                           _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                           _desc=domain( _name= ( boost::format( "%1%-%2%-%3%" ) % "hypercube" % convex_type().dimension() % 1 ).str() ,
                                         _shape="hypercube",
                                         _dim=convex_type().dimension(),
                                         _h=meshSize,
                                         _xmin=0,
                                         _xmax=L,
                                         _ymin=-r,
                                         _ymax=r  ) );

    mesh->addMarkerName( "Inlet", 1, 1 );
    mesh->addMarkerName( "Outlet", 3, 1 );
    mesh->addMarkerName( "wall1", 2, 1 );
    mesh->addMarkerName( "wall2", 4, 1 );
    std::cout << "number of elements of 2D: " << mesh->numElements() << "\n";
    //************************ Cylindre ****************************************
#elif (STOKESPRESSMESHTYPE == 2 ) 
    GeoTool::Node Centre(0,0,0);
    GeoTool::Node Rayon( r );
    GeoTool::Node Dir(1,0,0);
    GeoTool::Node Lg(L,0,0);
    GeoTool::Cylindre C( meshSize,"Cyl",Centre,Dir,Rayon,Lg);
    C.setMarker(_type="surface",_name="Inlet",_marker1=true);
    C.setMarker(_type="surface",_name="Outlet",_marker2=true);
    C.setMarker(_type="surface",_name="Wall",_marker3=true);
    C.setMarker(_type="volume",_name="OmegaFluide",_markerAll=true);

    mesh = C.createMesh(_mesh= new mesh_type,
                        _name="mymesh",
                        _partitions=Environment::worldComm().localSize(),
                        _worldcomm=Environment::worldComm() );

    //*************************************************************************
#elif (STOKESPRESSMESHTYPE == 3)
   //************************ Aorte ****************************************
    mesh = loadGMSHMesh( _mesh=new mesh_type,
                         _filename="/scratch/projet10/celine/aorta_M2.msh",
                         _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES,
                         _physical_are_elementary_regions=false,
                         _rebuild_partitions=true,
                         _partitions=Environment::worldComm().localSize(),
                         _worldcomm=Environment::worldComm() );

    std::cout << "number of elements of 3D: " << mesh->numElements() << "\n";

    //*********************************************************************

#endif

    Xh = space_type::New( mesh );
    std::cout << "\nnDof "  << Xh->nDof() << std::endl;
    P7=space_type_U::New(mesh);
    F = M_backend->newVector( Xh );
    D =  M_backend->newMatrix( Xh, Xh );
}
void
Stokes_Neumann_Neumann::run()
{
    this->init();

    auto U = Xh->element( "(u,p)" );
    auto V = Xh->element( "(u,q)" );
    auto u = U.element<0>( "u" );
    auto v = V.element<0>( "u" );
    auto p = U.element<1>( "p" );
    auto q = V.element<1>( "p" );
#if defined( FEELPP_USE_LM )
    auto lambda = U.element<2>();
    auto nu = V.element<2>();
#endif
    //# endmarker4 #

    Log() << "Data Summary:\n";
    Log() << "   hsize = " << meshSize << "\n";
    Log() << "  export = " << this->vm().count( "export" ) << "\n";
    Log() << "      mu = " << mu << "\n";
    Log() << " bccoeff = " << penalbc << "\n";




    //# marker5 #
    auto deft = gradt( u );
    auto def = grad( v );
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
    // //*********************  solution exacte (profile de poiseuille)2D *****************
    // exact solution known for hypercube (2D toy model)
    auto u_exact=vec((1-(Py()*Py())/(r*r))*(P_inlet-P_outlet)/(2*mu*L),cst(0.) );

    auto p_exact=(P_outlet-P_inlet)*Px()/L + P_inlet;

    auto f=vec(cst(0.) , cst(0.) );
    //*************************************************************************************


#elif  (STOKESPRESSMESHTYPE == 2)
    //*********************  solution exacte (profile de poiseuille)3D **********************
    // exact solution known for cylinder (3D toy model)
    auto u_exact=vec(  (P_inlet-P_outlet)*r*r*(1-(Py()*Py()+Pz()*Pz())/(r*r))/(4*L) , cst(0.) , cst(0.) );

    auto p_exact=(-Px()*(P_inlet-P_outlet)/L + P_inlet);

    auto f=vec(cst(0.) , cst(0.), cst(0.) );
    //*************************************************************************************

#elif (STOKESPRESSMESHTYPE == 3)
    // no exact solution known for the aorta (3D)
    // test case with f = 0
    auto f=vec(cst(0.) , cst(0.), cst(0.) );
#endif
    std::cout << "Treating the right hand side \n";

    // right hand side
    auto stokes_rhs = form1( _test=Xh, _vector=F );
    stokes_rhs += integrate( elements( mesh ),inner( f,id( v ) ) );
#if (STOKESPRESSMESHTYPE==1) || (STOKESPRESSMESHTYPE == 2)
    stokes_rhs += integrate( markedfaces( mesh,"Inlet" ),inner(-P_inlet*N(),id( v )) );
#elif (STOKESPRESSMESHTYPE == 3)
    stokes_rhs += integrate( markedfaces( mesh,"inlet" ),inner(-P_inlet*N(),id( v )) );
#endif
    /*
     * Construction of the left hand side
     */
    //# marker7 #
    std::cout << "Treating the left  hand side \n";
    auto stokes = form2( _test=Xh, _trial=Xh, _matrix=D );
    boost::timer chrono;
    stokes += integrate( elements( mesh ), mu*inner( deft,def ) );


    std::cout << "mu*inner(deft,def): " << chrono.elapsed() << "\n";
    chrono.restart();
    stokes +=integrate( elements( mesh ), - div( v )*idt( p ) + divt( u )*id( q ) );
    //***cylindre
#if (STOKESPRESSMESHTYPE == 2)
    stokes +=integrate( markedfaces( mesh, "Wall" ), -inner( SigmaNt,id( v ) ) );
    stokes +=integrate( markedfaces( mesh, "Wall" ), -inner( SigmaN,id( u ) ) );
    stokes +=integrate( markedfaces( mesh, "Wall" ), +penalbc*inner( idt( u ),id( v ) )/hFace() );
#elif (STOKESPRESSMESHTYPE == 1)
    //***Hypercube
    stokes +=integrate( markedfaces( mesh, "wall1" ), -inner( SigmaNt,id( v ) ) );
    stokes +=integrate( markedfaces( mesh, "wall2" ), -inner( SigmaNt,id( v ) ) );
    stokes +=integrate( markedfaces( mesh, "wall1" ), -inner( SigmaN,idt( u ) ) );
    stokes +=integrate( markedfaces( mesh, "wall2" ), -inner( SigmaN,idt( u ) ) );
    stokes +=integrate( markedfaces( mesh, "wall1" ), +penalbc*inner( idt( u ),id( v ) )/hFace() );
    stokes +=integrate( markedfaces( mesh, "wall2" ), +penalbc*inner( idt( u ),id( v ) )/hFace() );
#elif (STOKESPRESSMESHTYPE == 3)
    stokes +=integrate( markedfaces( mesh, "wall" ), -inner( SigmaNt,id( v ) ) );
    stokes +=integrate( markedfaces( mesh, "wall" ), -inner( SigmaN,idt( u ) ) );
    stokes +=integrate( markedfaces( mesh, "wall" ), +penalbc*inner( idt( u ),id( v ) )/hFace() );
#endif
    std::cout << "(u,p): " << chrono.elapsed() << "\n";
    chrono.restart();


    //# endmarker7 #

    M_backend->solve( _matrix=D, _solution=U, _rhs=F );

#if 0
    U.save( _path="." );
    u.save( _path="." );
    p.save( _path="." );
    V.load( _path="." );
    v.load( _path="." );
    q.load( _path="." );
    std::cout << "||u-v||=" << ( u-v ).l2Norm() << "\n";
    std::cout << "||p-q||=" << ( p-q ).l2Norm() << "\n";
#endif
    std::cout << "Export results..." << std::endl;
#if (STOKESPRESSMESHTYPE == 1) || (STOKESPRESSMESHTYPE == 2)
    this->exportResults( u_exact, p_exact, U, V );
#elif (STOKESPRESSMESHTYPE == 3)
    this->exportResults( U, V );
#endif
    Log() << "[dof]         number of dof: " << Xh->nDof() << "\n";
    Log() << "[dof]    number of dof/proc: " << Xh->nLocalDof() << "\n";
    Log() << "[dof]      number of dof(U): " << Xh->functionSpace<0>()->nDof()  << "\n";
    Log() << "[dof] number of dof/proc(U): " << Xh->functionSpace<0>()->nLocalDof()  << "\n";
    Log() << "[dof]      number of dof(P): " << Xh->functionSpace<1>()->nDof()  << "\n";
    Log() << "[dof] number of dof/proc(P): " << Xh->functionSpace<1>()->nLocalDof()  << "\n";
} // Stokes::run


#if (STOKESPRESSMESHTYPE == 1) || (STOKESPRESSMESHTYPE == 2)
template<typename ExprUExact, typename ExprPExact>
void
    Stokes_Neumann_Neumann::exportResults( ExprUExact u_exact, ExprPExact p_exact,
                       element_type& U, element_type& V )
#elif (STOKESPRESSMESHTYPE == 3 )
void
    Stokes_Neumann_Neumann::exportResults(element_type& U, element_type& V )
#endif
{
    auto u = U.element<0>();
    auto p = U.element<1>();

    auto v = V.element<0>();
    auto q = V.element<1>();
#if defined( FEELPP_USE_LM )
    auto lambda = U.element<2>();
    auto nu = V.element<2>();
    Log() << "value of the Lagrange multiplier lambda= " << lambda( 0 ) << "\n";
    std::cout << "value of the Lagrange multiplier lambda= " << lambda( 0 ) << "\n";

#endif
#if (STOKESPRESSMESHTYPE == 1) || (STOKESPRESSMESHTYPE == 2)
    double outflow_inlet = integrate( markedfaces( u.mesh(),"Inlet" ), inner( idv(u),N() ) ).evaluate()(0,0);
    double outflow_outlet = integrate( markedfaces( u.mesh(),"Outlet" ), inner( idv(u),N() ) ).evaluate()(0,0);
#elif (STOKESPRESSMESHTYPE == 3)
    double outflow_inlet = integrate( markedfaces( u.mesh(),"inlet" ), inner( idv(u),N() ) ).evaluate()(0,0);
    double outflow_outlet = integrate( markedfaces( u.mesh(),"outlets"), inner( idv(u),N() ) ).evaluate()(0,0);
#endif
    std::cout<<"outflow inlet = "<< outflow_inlet<<"\n";
    std::cout<<"outflow outlet = "<< outflow_outlet<<"\n";
    std::cout<<"flow difference = "<<outflow_inlet+outflow_outlet<<"\n";
   // Cylindre //
#if (STOKESPRESSMESHTYPE == 2)
    // exact outflow known in cylinder toy model (3D)
    auto r=1;
    auto L=5;
    auto P_inlet=1.;
    auto P_outlet=0.;
    double outflow_exact = ((4*math::atan(1.) * r * r * r * r)/(8*L))*(P_inlet-P_outlet);
    std::cout<<"exact outflow = "<<outflow_exact<<"\n";
#endif

#if (STOKESPRESSMESHTYPE==1) || (STOKESPRESSMESHTYPE==2)
    v = vf::project( u.functionSpace(), elements( u.mesh() ), u_exact );
    q = vf::project( p.functionSpace(), elements( p.mesh() ), p_exact );

#endif
    if ( exporter->doExport() )
    {
        exporter->step( 0 )->setMesh( U.functionSpace()->mesh() );
        exporter->step( 0 )->addRegions();
        //auto v = U.functionSpace()->functionSpace<0> ()->element();
        //v = U.element<0>();
        exporter->step( 0 )->add( "u", U.element<0>() );
        exporter->step( 0 )->add( "p", U.element<1>() );
#if (STOKESPRESSMESHTYPE==1)|| (STOKESPRESSMESHTYPE == 2)
        exporter->step( 0 )->add( "u_exact", V.element<0>() );
        exporter->step( 0 )->add( "p_exact", V.element<1>() );
#endif
        exporter->save();
    }

} // Stokes::export
} // Feel

int
main( int argc, char** argv )
{

    using namespace Feel;

    Environment env( argc, argv );

    /* assertions handling */
    Feel::Assert::setLog( "stokes.assert" );

    /* define and run application */
    Feel::Stokes_Neumann_Neumann Stokes_Neumann_Neumann( argc, argv, makeAbout(), makeOptions() );
    Stokes_Neumann_Neumann.run();
}
