/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   Date: 2005-11-08

   Copyright (C) 2005,2006 EPFL
   Copyright (C) 2007 Université Joseph Fourier (Grenoble I)

   This library is free softwarey; you can redistribute it and/or
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
   \file ale.hpp
   \author Goncalo Pena <goncalo.pena@epfl.ch>
   \date 2007-12-12
*/



#ifndef __TestALE
#define __TestALE 1

#include <feel/options.hpp>

#include <feel/feelpoly/fekete.hpp>
#include <feel/feelalg/backend.hpp>

#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelfilters/gmshhypercubedomain.hpp>

#include <feel/feelcore/applicationxml.hpp>

#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>

#include <feel/feeldiscr/operatorlagrangep1.hpp>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/interpolate.hpp>


#include <feel/feeldiscr/ale.hpp>


inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description meshHighOrderoptions( "TestALE options" );
    meshHighOrderoptions.add_options()
    ( "hsize", Feel::po::value<double>()->default_value( 2 ), "meshsize" )
    ;

    return meshHighOrderoptions.add( Feel::feel_options() );
}

inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "TestALE" ,
                           "TestALE" ,
                           "0.2",
                           "Test High order meshes",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2007 Universite Joseph Fourier" );

    about.addAuthor( "Goncalo Pena", "developer", "gpena@mat.uc.pt", "" );
    return about;
}


namespace Feel
{
template< int N >
class TestALE
    :
public ApplicationXML
{
    typedef ApplicationXML super;

    typedef Simplex<2, 1> convex_type;

    typedef Mesh<Simplex<1, 1> > struct_mesh_type;
    typedef Mesh<Simplex<2, 1> > mesh_type;
    typedef Mesh<Simplex<2, N> > new_mesh_type;

    typedef bases<Lagrange<N, Scalar, Continuous, PointSetFekete> > struct_basis_type;
    typedef FunctionSpace<struct_mesh_type, struct_basis_type, double> struct_functionspace_type;
    typedef boost::shared_ptr<struct_functionspace_type> struct_functionspace_ptrtype;
    typedef typename struct_functionspace_type::element_type struct_element_type;

    typedef typename PointSetEquiSpaced<Hypercube<1,1>, N, double>::points_type node_points_type;

    typedef bases<Lagrange<1, Vectorial, Continuous, PointSetFekete> > p1_ale_basis_type;
    typedef FunctionSpace< mesh_type, p1_ale_basis_type, double> p1_functionspace_type;
    typedef boost::shared_ptr<p1_functionspace_type> p1_functionspace_ptrtype;
    typedef typename p1_functionspace_type::element_type p1_element_type;

    typedef bases<Lagrange<N, Vectorial, Continuous, PointSetFekete> > pN_ale_basis_type;
    typedef FunctionSpace< new_mesh_type, pN_ale_basis_type, double> pN_visualize_functionspace_type;
    typedef boost::shared_ptr<pN_visualize_functionspace_type> pN_visualize_functionspace_ptrtype;
    typedef typename pN_visualize_functionspace_type::element_type pN_element_type;


    typedef bases<Lagrange<N, Scalar, Continuous, PointSetFekete> > basis_type;
    typedef FunctionSpace< new_mesh_type, basis_type, double> functionspace_type;
    typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;
    typedef typename functionspace_type::element_type fs_element_type;


    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;

    /*matrix*/
    typedef Backend<double> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

public:

    TestALE( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) ),
        M_backend( backend_type::build( this->vm() ) )
    {
        meshSize = this->vm()["hsize"].template as<double>();

        Parameter h;

        if ( N < 4 )
            h=Parameter( _name="h",_type=CONT_ATTR,_cmdName="hsize",_values="0.1:0.05:0.5" );

        else
            h=Parameter( _name="h",_type=CONT_ATTR,_cmdName="hsize",_values="0.5:0.5:1" );

        this->addParameter( Parameter( _name="order",_type=DISC_ATTR,_values=boost::lexical_cast<std::string>( N  ).c_str() ) )
        .addParameter( h );

        std::vector<Parameter> depend;
        std::vector<std::string> funcs;
        depend.push_back( h );
        std::ostringstream oss;
        oss << "h**" << boost::lexical_cast<std::string>( N + 1  ) ;
        funcs.push_back( oss.str() );

        this->
        addOutput( Output( _name="norm_L2",_latex="\\left\\| . \\right\\|_{L^2}",_dependencies=depend,_funcs=funcs ) );
        std::cout << "do_export = " << exporter->doExport() << "\n";
    }

    /**
     * alias for run()
     */
    void operator()()
    {
        run();
    }

    void exportResults( double tn, pN_element_type& U )
    {
        if ( exporter->doExport() )
        {
#if 1
            functionspace_ptrtype Wh = functionspace_type::New( U.functionSpace()->mesh() );

            OperatorLagrangeP1<functionspace_type> I( Wh, M_backend );

            typename OperatorLagrangeP1<functionspace_type>::dual_image_space_ptrtype Yh( I.dualImageSpace() );
            exporter->step( tn )->setMesh( Yh->mesh() );

#else
            exporter->step( tn )->setMesh( U.mesh() );
#endif
            exporter->step( tn )->add( "U", U );
            exporter->save();
        }
    }

    /**
     * run the convergence test
     */
    void run( );

private:

    double meshSize;

    export_ptrtype exporter;

    backend_ptrtype M_backend;
}; // TestALE



template<int N>
void
TestALE<N>::run()
{
    this->addParameterValue( N )
    .addParameterValue( this->vm()["hsize"].template as<double>() );

    if ( this->preProcessing() == RUN_EXIT ) return;

    boost::timer time;
    using namespace Feel::vf;


    /*
      Define mesh for the curved boundaries
    */
    boost::shared_ptr<struct_mesh_type> struct_mesh( new struct_mesh_type );
    GmshHypercubeDomain td( 1,1,1,false );
    td.setCharacteristicLength( meshSize );
    td.setX( std::make_pair( 0, 5 ) );
    std::string fname = td.generate( "boundary_description" );
    ImporterGmsh<struct_mesh_type> struct_import( fname );
    struct_mesh->accept( struct_import );

    std::cout << "mesh generation in " << time.elapsed() << "s\n";
    time.restart();

    /*
      Define functionspace for the boundaries as well as the polynomials that define it
    */
    struct_functionspace_ptrtype Xh = struct_functionspace_type::New( struct_mesh );
    struct_element_type bc_top( Xh, "bc_top" );
    AUTO( f,  val( 1.0 + 0.3*cos( Px() ) )  );
    bc_top = vf::project( Xh, elements( Xh->mesh() ), f );

    std::cout << "bc top in " << time.elapsed() << "s\n";
    time.restart();

    struct_element_type bc_bottom( Xh, "bc_bottom" );

    AUTO( f2, val( -1.1 - 0.3*cos( Px() ) )  );

    bc_bottom = vf::project( Xh, elements( Xh->mesh() ), f2 );

    std::cout << "bc bottom in " << time.elapsed() << "s\n";
    time.restart();

    struct_element_type bc_reference_top( Xh, "bc_reference_top" );
    bc_reference_top = vf::project( Xh, elements( Xh->mesh() ), constant( 1.0 ) );

    std::cout << "ref bc top in " << time.elapsed() << "s\n";
    time.restart();

    struct_element_type bc_reference_bottom( Xh, "bc_reference_bottom" );
    bc_reference_bottom = vf::project( Xh, elements( Xh->mesh() ), constant( -1.0 ) );

    std::cout << "ref bc bottom in " << time.elapsed() << "s\n";
    time.restart();

    struct_element_type bc_disp_top( Xh, "bc_top" );
    bc_disp_top = bc_top;
    bc_disp_top -= bc_reference_top;
    bc_disp_top.updateGlobalValues();

    struct_element_type bc_disp_bottom( Xh, "bc_bottom" );
    bc_disp_bottom = bc_bottom;
    bc_disp_bottom -= bc_reference_bottom;
    bc_disp_bottom.updateGlobalValues();

    std::cout << "bc disp top/bottom in " << time.elapsed() << "s\n";
    time.restart();
    /*
      Define mesh for the domain
    */
    boost::shared_ptr<mesh_type> mesh( new mesh_type );

    GmshHypercubeDomain td2( 2,1,2,false );
    td2.setCharacteristicLength( meshSize );
    td2.setX( std::make_pair( 0, 5 ) );
    td2.setY( std::make_pair( -1, 1 ) );
    std::string fname2 = td2.generate( "geometry" );
    ImporterGmsh<mesh_type> import( fname2 );
    mesh->accept( import );

    std::cout << "Ah mesh in " << time.elapsed() << "s\n";
    time.restart();

    p1_functionspace_ptrtype Ah = p1_functionspace_type::New( mesh );


    //define set of flags
    std::map< std::string, std::vector<flag_type> > flagSet;
    flagSet["fixed_bc"].push_back( 3 );
    flagSet["fixed_bc"].push_back( 1 );
    flagSet["moving_bc"].push_back( 2 );
    flagSet["moving_bc"].push_back( 4 );

    //define set of polynomials that describe boundary
    std::vector<struct_element_type> polyBoundarySet;
    polyBoundarySet.push_back( bc_disp_bottom );
    polyBoundarySet.push_back( bc_disp_top );


    std::vector<struct_element_type> referencePolyBoundarySet;
    referencePolyBoundarySet.push_back( bc_reference_bottom );
    referencePolyBoundarySet.push_back( bc_reference_top );

    std::cout << "Ah space in " << time.elapsed() << "s\n";
    time.restart();

    /*
      Create the ale map
    */
    ALE< Simplex<2,N> > aleFactory( std::make_pair( 0,5 ), mesh, this->vm() );

    aleFactory.generateHighOrderMap( flagSet["moving_bc"], referencePolyBoundarySet, polyBoundarySet );

    std::cout << "ALE ho map in " << time.elapsed() << "s\n";
    time.restart();

    MeshHighOrder< Simplex<2, N> > auxiliar_mesh ( mesh );
    auxiliar_mesh.generateMesh( flagSet["moving_bc"], referencePolyBoundarySet );
    boost::shared_ptr<new_mesh_type> aux_mesh = auxiliar_mesh.getMesh();

    std::cout << "ALE ho map mesh in " << time.elapsed() << "s\n";
    time.restart();

    pN_visualize_functionspace_ptrtype visH = pN_visualize_functionspace_type::New( aux_mesh );
    pN_element_type aux_element( visH, "aux" );
    aux_element = vf::project( visH, elements( visH->mesh() ), sin( Px() )*cos( Py() )*oneX() + sin( Py() )*cos( Px() )*oneY() );
    this->exportResults( 0.0, aux_element );

    std::cout << "ALE visu in " << time.elapsed() << "s\n";
    time.restart();
    /*
      Test ale map in boundary
    */
    double error_bottom_first = math::sqrt( integrate( markedfaces( Ah->mesh(), 2 ),
                                            ( trans( idv( aleFactory.getMap() ) )*oneX() - Px() )*( trans( idv( aleFactory.getMap() ) )*oneX() - Px() ),
                                            _Q<10>()
                                                     ).evaluate()( 0, 0 ) );

    double error_bottom_second = math::sqrt( integrate( markedfaces( Ah->mesh(), 2 ),
                                 ( trans( idv( aleFactory.getMap() ) )*oneY() - f2 )*( trans( idv( aleFactory.getMap() ) )*oneY() - f2 ),
                                 _Q<10>()
                                                      ).evaluate()( 0, 0 ) );

    std::cout << "error bottom in " << time.elapsed() << "s\n";
    time.restart();
    std::cout << "Error in first component of ALE map: " << error_bottom_first << "\n";
    std::cout << "Error in second component of ALE map: " << error_bottom_second << "\n";



    double error_bottom = math::sqrt( integrate( markedfaces( Ah->mesh(), 2 ),
                                      trans( idv( aleFactory.getMap() ) - vec( Px(),f2 ) )*( idv( aleFactory.getMap() ) - vec( Px(),f2 ) ),
                                      _Q<10>()
                                               ).evaluate()( 0, 0 ) );

    double error_top = math::sqrt( integrate( markedfaces( Ah->mesh(), 4 ),
                                   trans( idv( aleFactory.getMap() ) - vec( Px(),f ) )*( idv( aleFactory.getMap() ) - vec( Px(),f ) ),
                                   _Q<10>()
                                            ).evaluate()( 0, 0 ) );

    std::cout << "error top in " << time.elapsed() << "s\n";
    time.restart();
    std::cout << "Error top in the boundary: " << error_top << "\n";
    std::cout << "Error bottom in the boundary: " << error_bottom << "\n";
    double errbdy = math::sqrt( error_top*error_top + error_bottom*error_bottom );
    std::cout << "Error in the boundary: " <<  errbdy << "\n";
    Log() << "Error in the boundary: " << errbdy << "\n";


    MeshMover<new_mesh_type> mesh_mover;
    mesh_mover.apply( visH->mesh(), aleFactory.getDisplacement() );

    std::cout << "mesh move in " << time.elapsed() << "s\n";
    time.restart();

    aux_element = vf::project( visH, elements( visH->mesh() ), sin( Px() )*cos( Py() )*oneX() + sin( Py() )*cos( Px() )*oneY() );
    aux_element.updateGlobalValues();
    this->exportResults( 1.0, aux_element );

    std::cout << "mesh moved visu in " << time.elapsed() << "s\n";
    time.restart();

    this->addOutputValue( errbdy );
    this->postProcessing();

} // end run routine

} // end Feel
#endif // __TestALE
