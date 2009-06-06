/* -*- mode: c++ -*-

   This file is part of the Life library

   Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   Date: 2005-11-08

   Copyright (C) 2005,2006 EPFL
   Copyright (C) 2007 Université Joseph Fourier (Grenoble I)

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
   \file test_mesh.hpp
   \author Goncalo Pena <goncalo.pena@epfl.ch>
   \date 2007-12-12
*/



#ifndef __TestALE
#define __TestALE 1

#include <life/options.hpp>

#include <life/lifepoly/fekete.hpp>
#include <life/lifealg/backendtrilinos.hpp>
#include <life/lifealg/backend.hpp>

#include <life/lifefilters/importergmsh.hpp>
#include <life/lifefilters/gmsh.hpp>
#include <life/lifefilters/gmshtensorizeddomain.hpp>

#include <life/lifecore/application.hpp>

#include <life/lifevf/vf.hpp>
#include <life/lifefilters/exporter.hpp>
#include <life/lifefilters/exporterensight.hpp>

#include <life/lifediscr/operatorlagrangep1.hpp>

#include <life/lifediscr/mesh.hpp>
#include <life/lifediscr/functionspace.hpp>
#include <life/lifediscr/interpolate.hpp>


#include <life/lifediscr/ale.hpp>


inline
Life::po::options_description
makeOptions()
{
    Life::po::options_description meshHighOrderoptions("TestALE options");
    meshHighOrderoptions.add_options()
        ("h", Life::po::value<double>()->default_value( 2 ), "meshsize")
        ;

    return meshHighOrderoptions.add( Life::life_options() );
}

inline
Life::AboutData
makeAbout()
{
    Life::AboutData about( "TestALE" ,
                           "TestALE" ,
                           "0.1",
                           "Test High order meshes",
                           Life::AboutData::License_GPL,
                           "Copyright (c) 2007 Universite Joseph Fourier");

    about.addAuthor("Goncalo Pena", "developer", "goncalo.pena@epfl.ch", "");
    return about;
}


namespace Life
{
template< int N >
class TestALE
    :
        public Application
{
    typedef Application super;

    typedef Simplex<2, 1> convex_type;

    typedef Mesh<GeoEntity<Simplex<1, 1> > > struct_mesh_type;
    typedef Mesh<GeoEntity<Simplex<2, 1> > > mesh_type;
    typedef Mesh<GeoEntity<Simplex<2, N> > > new_mesh_type;

    typedef fusion::vector<fem::Lagrange<1, N, Scalar, Continuous, double, Simplex, PointSetWarpBlend> > struct_basis_type;
    typedef FunctionSpace<struct_mesh_type, struct_basis_type, double> struct_functionspace_type;
    typedef boost::shared_ptr<struct_functionspace_type> struct_functionspace_ptrtype;
    typedef typename struct_functionspace_type::element_type struct_element_type;

    typedef typename PointSetEquiSpaced<SimplexProduct<1,1>, N, double>::points_type node_points_type;

    typedef fusion::vector<fem::Lagrange<2, 1, Vectorial, Continuous, double, Simplex, PointSetFekete> > p1_ale_basis_type;
    typedef FunctionSpace< mesh_type, p1_ale_basis_type, double> p1_functionspace_type;
    typedef boost::shared_ptr<p1_functionspace_type> p1_functionspace_ptrtype;

    typedef fusion::vector<fem::Lagrange<2, N, Vectorial, Continuous, double, Simplex, PointSetFekete> > pN_ale_basis_type;
    typedef FunctionSpace< new_mesh_type, pN_ale_basis_type, double> pN_visualize_functionspace_type;
    typedef boost::shared_ptr<pN_visualize_functionspace_type> pN_visualize_functionspace_ptrtype;
    typedef typename pN_visualize_functionspace_type::element_type pN_element_type;


    typedef fusion::vector<fem::Lagrange<2, N, Scalar, Continuous, double, Simplex, PointSetFekete> > basis_type;
    typedef FunctionSpace< new_mesh_type, basis_type, double> functionspace_type;
    typedef boost::shared_ptr<functionspace_type> functionspace_ptrtype;
    typedef typename functionspace_type::element_type fs_element_type;


    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;
    typedef typename export_type::timeset_type timeset_type;




    /*matrix*/
    typedef Backend<double> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    typedef OperatorMatrix op_mat_type;
    typedef boost::shared_ptr<OperatorMatrix> op_mat_ptrtype;

    typedef IM<2, 3*N, double, Simplex> im_type;





public:

    TestALE( int argc, char** argv, AboutData const& ad )
        :
        super( argc, argv, ad ),
        exporter( new ExporterEnsight<mesh_type>( "testALE" ) ),
        timeSet( new timeset_type( "testALE" ) ),
        M_backend( backend_type::build( this->vm() ) )
    {
        meshSize = this->vm()["h"].template as<double>();
        timeSet->setTimeIncrement( 1.0 );
        exporter->addTimeSet( timeSet );
    }

    TestALE( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        exporter( new ExporterEnsight<mesh_type>( "testALE" ) ),
        timeSet( new timeset_type( "testALE" ) ),
        M_backend( backend_type::build( this->vm() ) )
    {
        meshSize = this->vm()["h"].template as<double>();
        timeSet->setTimeIncrement( 1.0 );
        exporter->addTimeSet( timeSet );
    }

    TestALE( TestALE const& tc )
        :
        super( tc ),
        meshSize( tc.meshSize ),
        exporter( new ExporterEnsight<mesh_type>( "testALE" ) ),
        timeSet( new timeset_type( "testALE" ) ),
        M_backend( tc.M_backend )
    {
        timeSet->setTimeIncrement( 1.0 );
        exporter->addTimeSet( timeSet );
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
        typename timeset_type::step_ptrtype timeStep = timeSet->step( tn );

        functionspace_ptrtype Wh = functionspace_type::New( U.functionSpace()->mesh() );

        OperatorLagrangeP1<functionspace_type> I( Wh, M_backend );

        typename OperatorLagrangeP1<functionspace_type>::dual_image_space_ptrtype Yh( I.dualImageSpace() );
        typename OperatorLagrangeP1<functionspace_type>::dual_image_space_type::element_type w( Yh, "w" );
        typename OperatorLagrangeP1<functionspace_type>::dual_image_space_type::element_type e( Yh, "e" );

        timeStep->setMesh( Yh->mesh() );

        timeStep->add( "compX", U.comp(X) );
        timeStep->add( "compY", U.comp(Y) );

        exporter->save();
    }

    /**
     * run the convergence test
     */
    void run( );

private:

    double meshSize;

    export_ptrtype exporter;
    typename export_type::timeset_ptrtype timeSet;

    backend_ptrtype M_backend;
}; // TestALE



template<int N>
void
TestALE<N>::run()
{
    this->changeRepository( boost::format( "%1%/P%2%/h_%3%/" )
                            % this->about().appName()
                            % N
                            % this->vm()["h"].template as<double>()
                            );

    using namespace Life::vf;

    boost::timer time;

    /*
      Define mesh for the curved boundaries
    */
    boost::shared_ptr<struct_mesh_type> struct_mesh( new struct_mesh_type );
    std::string fname;
    GmshTensorizedDomain<1,1,Simplex> td;
    td.setCharacteristicLength( meshSize );
    td.setX( std::make_pair( 0, 5 ) );
    fname = td.generate( "structure" );
    ImporterGmsh<struct_mesh_type> struct_import( fname );
    struct_mesh->accept( struct_import );


    /*
      Define functionspace for the boundaries as well as the polynomials that define it
    */
    struct_functionspace_ptrtype Xh = struct_functionspace_type::New( struct_mesh );
    struct_element_type bc_top( Xh, "bc_top" );
    AUTO( f,  (0.02096834639998*Px()*Px()*Px()    -0.09351398506873*Px()*Px()    -0.09961900352798*Px()+    1.3000 )  ) ;
    bc_top = project( Xh, elements(*struct_mesh), f );

    struct_element_type bc_bottom( Xh, "bc_bottom" );
    AUTO( f2, (-0.0210*Px()*Px()*Px() + 0.0935*Px()*Px() + 0.0996*Px() -1.4000)  );
    bc_bottom = project( Xh, elements(*struct_mesh), f2 );

    struct_element_type bc_reference_top( Xh, "bc_reference_top" );
    bc_reference_top = project( Xh, elements(*struct_mesh), constant(1.0) );

    struct_element_type bc_reference_bottom( Xh, "bc_reference_bottom" );
    bc_reference_bottom = project( Xh, elements(*struct_mesh), constant(-1.0) );

    struct_element_type bc_disp_top( Xh, "bc_top" );
    bc_disp_top = bc_top;
    bc_disp_top -= bc_reference_top;

    struct_element_type bc_disp_bottom( Xh, "bc_bottom" );
    bc_disp_bottom = bc_bottom;
    bc_disp_bottom -= bc_reference_bottom;



    /*
      Define mesh for the domain
    */
    boost::shared_ptr<mesh_type> mesh( new mesh_type );
    std::ostringstream ostr;

    ostr << "Mesh.MshFileVersion = 1;\n"
         << "h=" << meshSize << ";\n"
         << "Point(1) = { 0, -1,0.0,h};\n"
         << "Point(2) = { 5, -1,0.0,h};\n"
         << "Point(3) = { 5,  1,0.0,h};\n"
         << "Point(4) = { 0,  1,0.0,h};\n"
         << "Line(1) = {1,2};\n"
         << "Line(2) = {2,3};\n"
         << "Line(3) = {3,4};\n"
         << "Line(4) = {4,1};\n"
         << "Line Loop(7) = {4,3,2,1};\n"
         << "Plane Surface(8) = {7};\n"
         << "Physical Line(1) = {1};\n"
         << "Physical Line(2) = {2};\n"
         << "Physical Line(3) = {3};\n"
         << "Physical Line(4) = {4};\n"
         << "Physical Surface(8) = {8};\n";

    Gmsh gmsh;

    std::string meshName = gmsh.generate( "geometry", ostr.str() );

    time.restart();
    ImporterGmsh<mesh_type> import( meshName );

    mesh->accept( import );
    mesh->components().set( MESH_CHECK | MESH_RENUMBER | MESH_UPDATE_EDGES | MESH_UPDATE_FACES );
    mesh->updateForUse();





    p1_functionspace_ptrtype Ah = p1_functionspace_type::New( mesh );

    //define set of flags
    std::map< std::string, std::vector<flag_type> > flagSet;
    flagSet["fixed_bc"].push_back(2);
    flagSet["fixed_bc"].push_back(4);
    flagSet["moving_bc"].push_back(1);
    flagSet["moving_bc"].push_back(3);

    //define set of polynomials that describe boundary
    std::vector<struct_element_type> polyBoundarySet;
    polyBoundarySet.push_back(bc_disp_bottom);
    polyBoundarySet.push_back(bc_disp_top);


    std::vector<struct_element_type> referencePolyBoundarySet;
    referencePolyBoundarySet.push_back(bc_reference_bottom);
    referencePolyBoundarySet.push_back(bc_reference_top);


    /*
      Create the ale map
    */
    ALE< Simplex<2,N> > aleFactory( std::make_pair(0,5), mesh );
    aleFactory.generateHighOrderMap( flagSet["moving_bc"], referencePolyBoundarySet, polyBoundarySet );

    MeshHighOrder< Simplex<2, N> > auxiliar_mesh ( mesh );
    auxiliar_mesh.generateMesh(flagSet["moving_bc"], referencePolyBoundarySet);
    boost::shared_ptr<new_mesh_type> aux_mesh = auxiliar_mesh.getMesh();


    pN_visualize_functionspace_ptrtype visH = pN_visualize_functionspace_type::New( aux_mesh );
    pN_element_type aux_element( visH, "aux");
    aux_element = project( visH, elements(visH->mesh()), sin(Px())*cos(Py())*oneX() + sin(Py())*cos(Px())*oneY());
    this->exportResults( 0.0, aux_element);



    /*
      Test ale map in boundary
    */

    im_type im;
    double error_bottom = math::sqrt(integrate( markedfaces( mesh, 1  ),
                                                im,
                                                trans(idv(aleFactory.getMap()) - (f2*oneY() + Px()*oneX()))*(idv(aleFactory.getMap()) - (f2*oneY() + Px()*oneX()))
                                                ).evaluate()( 0, 0 ));

    double error_top = math::sqrt(integrate( markedfaces( mesh, 3  ),
                                             im,
                                             trans(idv(aleFactory.getMap()) - (f*oneY() + Px()*oneX()))*(idv(aleFactory.getMap()) - (f*oneY() + Px()*oneX()))
                                             ).evaluate()( 0, 0 ));

    std::cout << "Error in the boundary: " << error_top + error_bottom << "\n";



    MeshMover<new_mesh_type> mesh_mover;
    mesh_mover.apply(aux_mesh, aleFactory.getDisplacement() );

    aux_element = project( visH, elements(visH->mesh()), sin(Px())*cos(Py())*oneX() + sin(Py())*cos(Px())*oneY());
    this->exportResults( 1.0, aux_element);

} // end run routine

} // end Life
#endif // __TestALE
