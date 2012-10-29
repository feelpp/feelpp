/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil;
c-basic-offset: 4; show-trailing-whitespace: t -*-
vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

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
   \file stokes_laplacian_inlet.cpp
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
#include <feel/feeldiscr/operatorinterpolation.hpp>
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
    Feel::po::options_description  stokesLaplacianInletOptions( "Stokes_laplacian_inlet options" );
    stokesLaplacianInletOptions.add_options()
        ( "penal", Feel::po::value<double>()->default_value( 0.5 ), "penalisation parameter" )
        ( "f", Feel::po::value<double>()->default_value( 0 ), "forcing term" )
        ( "mu", Feel::po::value<double>()->default_value( 1.0 ), "reaction coefficient component" )
        ( "hsize", Feel::po::value<double>()->default_value( 0.1 ), "first h value to start convergence" )
        ( "bctype", Feel::po::value<int>()->default_value( 0 ), "0 = strong Dirichlet, 1 = weak Dirichlet" )
        ( "bccoeff", Feel::po::value<double>()->default_value( 100.0 ), "coeff for weak Dirichlet conditions" )
        ( "generateMesh",Feel::po::value<bool>()->default_value( false ), "True for generating submesh" )
         ( "export-matlab", "export matrix and vectors in matlab" )
    ;
    return  stokesLaplacianInletOptions.add( Feel::feel_options() ) ;
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
    Feel::AboutData about( "stokes_laplacian_inlet" ,
                           "stokes_laplacian_inlet" ,
                           "0.1",
                           "Stokes equation on simplices or simplex products",
                           Feel::AboutData::License_GPL, "Copyright (c) 2009-2012 Université de Grenoble 1(Joseph Fourier)" );

    about.addAuthor( "Christophe Prud'homme", "developer","christophe.prudhomme@ujf-grenoble.fr", "" );
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
class Stokes_laplacian_inlet
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
    typedef Simplex<3> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;


    //# marker1 #,
    typedef Lagrange<2, Vectorial> basis_u_type;
    typedef Lagrange<1, Scalar> basis_p_type;
    typedef Lagrange<0, Scalar> basis_l_type;

    typedef bases<basis_u_type,basis_p_type> basis_type;

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

    //****************************************************************
    // MeshInlet
    typedef Simplex<2,1,3> inlet_convex_type;
    typedef Mesh<inlet_convex_type> meshInlet_type;
    typedef boost::shared_ptr<meshInlet_type> meshInlet_ptrtype;

    // Functions Spaces
    typedef Lagrange<1,Scalar> basis_u_InletType;
    typedef bases<basis_u_InletType> basis_inlet_type;
    typedef FunctionSpace<meshInlet_type, basis_inlet_type> space_InletType_U;

    typedef boost::shared_ptr<space_InletType_U> space_ptrInletType_U;
    typedef FunctionSpace<meshInlet_type, basis_inlet_type> space_InletType;

    typedef boost::shared_ptr<space_InletType> space_ptrtype2;
    typedef space_InletType::element_type element_type2;
    typedef Exporter<meshInlet_type> exportInlet_type;
    //****************************************************************
    //Interpolation //

    typedef Lagrange<1,Scalar> basis_u_InterpType;
    typedef bases<basis_u_InterpType> basis_interp_type;
    typedef FunctionSpace<mesh_type, basis_interp_type> space_InterpType_U;
    typedef boost::shared_ptr<space_InterpType_U> space_ptrInterpType_U;
    typedef FunctionSpace<mesh_type, basis_interp_type> space_InterpType;
    typedef boost::shared_ptr<space_InterpType> space_ptrtype3;
    typedef space_InterpType::element_type element_type3;
    FEELPP_DONT_INLINE
    Stokes_laplacian_inlet( int argc, char** argv, AboutData const& ad, po::options_description const& od );
    void exportResults(element_type& u, element_type& v );

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


private:

    backend_ptrtype M_backend;
    backend_ptrtype M_backendInlet;
    double meshSize;

    double mu;
    double penalbc;
    bool generateMesh;
    mesh_ptrtype mesh;
    meshInlet_ptrtype meshInlet;
    space_ptrtype Xh;
    space_ptrtype2 Xh2;
    space_ptrtype3 Xh3;
    space_ptrtype_U P7;
    sparse_matrix_ptrtype M,D,D2,M2;
    vector_ptrtype F,F2;

    boost::shared_ptr<export_type> exporter;
    boost::shared_ptr<exportInlet_type> exporterInletProblem;

}; // Stokes


    Stokes_laplacian_inlet:: Stokes_laplacian_inlet( int argc, char** argv, AboutData const& ad,po::options_description const& od )
    :
    super( argc, argv, ad, od ),
    M_backend( backend_type::build( this->vm() ) ),
    M_backendInlet( backend_type::build( this->vm() ) ),
    meshSize( this->vm()["hsize"].as<double>() ),
    mu( this->vm()["mu"].as<value_type>() ),
    penalbc( this->vm()["bccoeff"].as<value_type>() ),
    generateMesh(this->vm()["generateMesh"].as<bool>()),
    exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) ),
    exporterInletProblem( Exporter<meshInlet_type>::New( this->vm(), "meshInlet"
) )
{

}

void
Stokes_laplacian_inlet::init()
{
    if ( this->vm().count( "help" ) )
        {
        std::cout << this->optionsDescription() << "\n";
        return;
    }
    if ( this->vm().count( "nochdir" ) == false )
        this->changeRepository( boost::format("doc/tutorial/%1%/%2%/Part%6%/P%3%P%4%/h_%5%/" )
                                % this->about().appName()
                                % convex_type::name()
                                % basis_u_type::nOrder % basis_p_type::nOrder
                                % this->vm()["hsize"].as<double>()
                                % Environment::numberOfProcessors() );

    //********************* MESH GENERATION **********************************

    //************************ Aorte ****************************************
    mesh = loadGMSHMesh( _mesh=new mesh_type,
                         _filename="/scratch/ranine/aorta_M2.msh",
                         _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES,
                         //_rebuild_partitions=true,
                         _partitions=Environment::worldComm().localSize(),
                         _worldcomm=Environment::worldComm() );
    std::cout << "number of elements of 3D: " << mesh->numElements() << "\n";
    std::cout<<"arrived 1\n";

    //******************** create submesh for inlet face ****************/
    if (generateMesh == true){
        meshInlet = createSubmesh(mesh,  markedfaces( mesh,"inlet" ));
        std::cout << "number of elements of 2D in meshInlet: " << meshInlet->numElements() << "\n";
        saveGMSHMesh(_mesh=meshInlet, _filename="InletFace.msh");
    }
    else{
        meshInlet = loadGMSHMesh( _mesh=new meshInlet_type,
                                  _filename="/scratch/ranine/InletFace.msh",
                                  _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES,
                                  _physical_are_elementary_regions=false,
                                  //_rebuild_partitions=true,
                                  _partitions=Environment::worldComm().localSize(),
                                  _worldcomm=Environment::worldComm() );
        std::cout << "number of elements of 2D: " << meshInlet->numElements() <<"\n";
    }
    std::cout<<"arrived 2\n";
    Xh2= space_InletType::New( meshInlet );
    Xh3 = space_InterpType::New( mesh );
    Xh = space_type::New( mesh );
    std::cout<<"espaces de functions ok"<<std::endl;
    std::cout<<"NDof Xh2 = "<<Xh2->nDof()<<std::endl;
    std::cout<<"NDof Xh = "<<Xh->nDof()<<std::endl;
    F2 = M_backendInlet->newVector( Xh2 );
    D2 =  M_backendInlet->newMatrix( Xh2, Xh2 );
    std::cout<<"F2,D2"<<std::endl;
    //  P7=space_type_U::New(mesh);
    F = M_backend->newVector( Xh );
    D =  M_backend->newMatrix( Xh, Xh );
    std::cout<<"F,D"<<std::endl;
}
void
Stokes_laplacian_inlet::run()
{
    this->init();
    std::cout<<"3 : fin init"<<std::endl;
    auto u2 = Xh2->element( "u" );
    auto u3 = Xh3->element();
    auto v2 = Xh2->element( "u" );
    auto f2 = cst(10.);
    std::cout<<"4: u2,u3,v2 "<<std::endl;
    //****************************************************************
    //inlet velocity : solving laplacian problem
    //****************************************************************

    auto b = form1( _test=Xh2, _vector=F2 );
    b = integrate( elements( meshInlet ),f2*id( v2 ) );
 
    std::cout<<"5 : assemble form1 => ok"<<std::endl;
    auto a = form2( _test=Xh2, _trial=Xh2, _matrix=D2 );
    a = integrate( elements( meshInlet ), gradt( u2 )*trans( grad( v2 ) ) );
    D2->close();
    a += on( _range=boundaryfaces( meshInlet ),
             _element=u2,_rhs=F2,_expr=cst(0.) );
    std::cout<<"6 : assemble form2 => ok"<<std::endl;

    M_backendInlet->solve( _matrix=D2, _solution=u2, _rhs=F2 );
    //****************************************************************
    exporterInletProblem->step( 0 )->setMesh( u2.functionSpace()->mesh() );
    exporterInletProblem->step( 0 )->add( "u2",u2);
    exporterInletProblem->save();
    std::cout<<"7 : exporter Probleme inlet => ok"<<std::endl;
    //****************************************************************
    // interpolation : imposer la vitesse en entrée
 
    auto opI = opInterpolation(_domainSpace = Xh2,
                               _imageSpace = Xh3,
                               _type = InterpolationConforme(),
                               _range =  markedfaces( mesh,"inlet" ));
    std::cout<<"8 : interpolation 1/2"<<std::endl;
    opI->apply(u2 ,u3);
    std::cout<<"9 : interpolation 2/2"<<std::endl;
    auto U = Xh->element( "(u,p)" );
    auto V = Xh->element( "(u,q)" );
    auto u = U.element<0>( "u" );
    auto v = V.element<0>( "u" );
    auto p = U.element<1>( "p" );
    auto q = V.element<1>( "p" );
    auto P_inlet=1;
    auto P_outlet=0.;
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
    auto f=vec(cst(0.) , cst(0.), cst(0.) );
    std::cout << "Treating the right hand side \n";
    // right hand side
    auto stokes_rhs = form1( _test=Xh, _vector=F );
    stokes_rhs += integrate( markedfaces( mesh,"inlet" ),inner(-abs(idv(u3)/cst(20.))*N() ,-SigmaN+penalbc*id( v )/hFace() ) );

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
    stokes +=integrate( elements( mesh ), - div( v )*idt( p ) + divt( u )*id( q) );
    stokes +=integrate( markedfaces( mesh, "wall" ), -inner( SigmaNt,id( v ) ));
    stokes +=integrate( markedfaces( mesh, "wall" ), -inner( SigmaN,idt( u ) ));
    stokes +=integrate( markedfaces( mesh, "wall" ), +penalbc*inner( idt( u),id( v ) )/hFace() );
    stokes +=integrate( markedfaces( mesh, "inlet" ), -inner( SigmaNt,id( v ) ));
    stokes +=integrate( markedfaces( mesh, "inlet" ), -inner( SigmaN,idt( u ) ));
    stokes +=integrate( markedfaces( mesh, "inlet" ), +penalbc*inner( idt( u),id( v ) )/hFace() );
    std::cout << "(u,p): " << chrono.elapsed() << "\n";
    chrono.restart();
    M_backend->solve( _matrix=D, _solution=U, _rhs=F );
    std::cout<<"resolution => ok"<<std::endl;
    this->exportResults( U, V );
}
void
Stokes_laplacian_inlet::exportResults(element_type& U, element_type& V ){
    auto u = U.element<0>();
    auto p = U.element<1>();

    auto v = V.element<0>();
    auto q = V.element<1>();
    double outflow_inlet = integrate( markedfaces( u.mesh(),"inlet" ), inner(idv(u),N() ) ).evaluate()(0,0);
    double outflow_outlet = integrate( markedfaces( u.mesh(),"outlets"), inner(idv(u),N() ) ).evaluate()(0,0);
    std::cout<<"outflow inlet = "<< outflow_inlet<<"\n";
    std::cout<<"outflow outlet = "<< outflow_outlet<<"\n";
    std::cout<<"flow difference = "<<outflow_inlet+outflow_outlet<<"\n";
    if ( exporter->doExport() )
    {
        exporter->step( 0 )->setMesh( U.functionSpace()->mesh() );
        exporter->step( 0 )->addRegions();
        exporter->step( 0 )->add( "u", U.element<0>() );
        exporter->step( 0 )->add( "p", U.element<1>() );
    }
        exporter->save();
}
}
int
main( int argc, char** argv )
{

    using namespace Feel;

    Environment env( argc, argv );

    /* assertions handling */
    Feel::Assert::setLog( "stokes.assert" );

    // SOME BAD ELEMENTS
    // P1/P0 : locking
    //typedef Feel::Stokes<Simplex<2>, Lagrange<1, Vectorial>,Lagrange<0,Scalar,Discontinuous> > stokes_type;
    // P1/P1 : spurious modes
    //typedef Feel::Stokes<Simplex<2>, Lagrange<1, Vectorial>,Lagrange<1,Scalar> > stokes_type;

    // SOME GOOD ELEMENTS
    // P2/P1
    // CR0/P0
    //typedef Feel::Stokes<Simplex<2>, CrouzeixRaviart<1, Vectorial>,Lagrange<0,Scalar,Discontinuous> > stokes_type;


    /* define and run application */
    Feel::Stokes_laplacian_inlet Stokes_laplacian_inlet( argc, argv, makeAbout(), makeOptions() );
    Stokes_laplacian_inlet.run();
}
