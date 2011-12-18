/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
  Author(s): Cecile Daversin  <cecile.daversin@lncmi.cnrs.fr>
       Date: 2011-12-07

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
   \file test_hcurl.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \author Cecile Daversin <cecile.daversin@lncmi.cnrs.fr>
   \date 2011-12-07
 */
#define USE_BOOST_TEST 1

// make sure that the init_unit_test function is defined by UTF
//#define BOOST_TEST_MAIN
// give a name to the testsuite
#define BOOST_TEST_MODULE H_curl approximation
// disable the main function creation, use our own
//#define BOOST_TEST_NO_MAIN

#if defined(USE_BOOST_TEST)
#include <boost/test/unit_test.hpp>
using boost::unit_test::test_suite;
#include <boost/test/floating_point_comparison.hpp>
#endif

#include <feel/feelcore/application.hpp>

/** include predefined feel command line options */
#include <feel/options.hpp>

/** include linear algebra backend */
#include <feel/feelalg/backend.hpp>

/** include linear algebra backend */
#include <feel/feelalg/vector.hpp>

/** include function space class */
#include <feel/feeldiscr/functionspace.hpp>

/** include helper function to define \f$P_0\f$ functions associated with regions  */
#include <feel/feeldiscr/region.hpp>

/** include integration methods */
#include <feel/feelpoly/im.hpp>

/** include gmsh mesh importer */
#include <feel/feelfilters/gmsh.hpp>

/** include exporter factory class */
#include <feel/feelfilters/exporter.hpp>

/** include  polynomialset header */
#include <feel/feelpoly/polynomialset.hpp>
#include <feel/feelpoly/lagrange.hpp>
#include <feel/feelpoly/nedelec.hpp>

/** include  the header for the variational formulation language (vf) aka FEEL++ */
#include <feel/feelvf/vf.hpp>

#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <feel/feelcore/traits.hpp>

#include <feel/feelalg/datamap.hpp>

namespace Feel
{
/// Geometry for one-element meshes
gmsh_ptrtype
oneelement_geometry_ref( double h = 2)
{
    std::ostringstream costr;
    costr <<"Mesh.MshFileVersion = 2.2;\n"
          <<"Mesh.CharacteristicLengthExtendFromBoundary=1;\n"
          <<"Mesh.CharacteristicLengthFromPoints=1;\n"
          <<"Mesh.ElementOrder=1;\n"
          <<"Mesh.SecondOrderIncomplete = 0;\n"
          <<"Mesh.Algorithm = 6;\n"
          <<"Mesh.OptimizeNetgen=1;\n"
          <<"// partitioning data\n"
          <<"Mesh.Partitioner=1;\n"
          <<"Mesh.NbPartitions=1;\n"
          <<"Mesh.MshFilePartitioned=0;\n"
          <<"h=" << h << ";\n"
          <<"Point(1) = {-1,-1,0,h};\n"
          <<"Point(2) = {1,-1,0,h};\n"
          <<"Point(3) = {-1,1,0,h};\n"
          <<"Line(1) = {1,2};\n"
          <<"Line(2) = {2,3};\n"
          <<"Line(3) = {3,1};\n"
          <<"Transfinite Line{1} = 1;\n"
          <<"Transfinite Line{2} = 1;\n"
          <<"Transfinite Line{3} = 1;\n"
          <<"Line Loop(4) = {3,1,2};\n"
          <<"Plane Surface(5) = {4};\n"
          <<"Physical Line(\"hor\") = {1};\n"
          <<"Physical Line(\"hypo\") = {2};\n"
          <<"Physical Line(\"vert\") = {3};\n"
          <<"Physical Surface(9) = {5};\n";

    std::ostringstream nameStr;
    if ( std::abs( h - 2 ) < 1e-10 )
        nameStr << "one-elt-ref";
    else
        nameStr << "one-elt-mesh-ref";
    gmsh_ptrtype gmshp( new Gmsh );
    gmshp->setPrefix( nameStr.str() );
    gmshp->setDescription( costr.str() );
    return gmshp;
}

// homothetic transformation of reference element (center 0, rate 2)
gmsh_ptrtype
oneelement_geometry_real_1(double h = 2)
{
    std::ostringstream costr;
    costr <<"Mesh.MshFileVersion = 2.2;\n"
          <<"Mesh.CharacteristicLengthExtendFromBoundary=1;\n"
          <<"Mesh.CharacteristicLengthFromPoints=1;\n"
          <<"Mesh.ElementOrder=1;\n"
          <<"Mesh.SecondOrderIncomplete = 0;\n"
          <<"Mesh.Algorithm = 6;\n"
          <<"Mesh.OptimizeNetgen=1;\n"
          <<"// partitioning data\n"
          <<"Mesh.Partitioner=1;\n"
          <<"Mesh.NbPartitions=1;\n"
          <<"Mesh.MshFilePartitioned=0;\n"
          <<"h=" << h << ";\n"
          <<"Point(1) = {-2,-2,0,h};\n"
          <<"Point(2) = {2,-2,0,h};\n"
          <<"Point(3) = {-2,2,0,h};\n"
          <<"Line(1) = {1,2};\n"
          <<"Line(2) = {2,3};\n"
          <<"Line(3) = {3,1};\n"
          <<"Transfinite Line{1} = 1;\n"
          <<"Transfinite Line{2} = 1;\n"
          <<"Transfinite Line{3} = 1;\n"
          <<"Line Loop(4) = {3,1,2};\n"
          <<"Plane Surface(5) = {4};\n"
          <<"Physical Line(\"hor\") = {1};\n"
          <<"Physical Line(\"hypo\") = {2};\n"
          <<"Physical Line(\"vert\") = {3};\n"
          <<"Physical Surface(9) = {5};\n";

    std::ostringstream nameStr;
    if ( std::abs( h - 2 ) < 1e-10 )
        nameStr << "one-elt-real-homo";
    else
        nameStr << "one-elt-mesh-homo";
    gmsh_ptrtype gmshp( new Gmsh );
    gmshp->setPrefix( nameStr.str() );
    gmshp->setDescription( costr.str() );
    return gmshp;
}

// Rotation of angle (pi/2)
gmsh_ptrtype
oneelement_geometry_real_2(double h = 2)
{
    std::ostringstream costr;
    costr <<"Mesh.MshFileVersion = 2.2;\n"
          <<"Mesh.CharacteristicLengthExtendFromBoundary=1;\n"
          <<"Mesh.CharacteristicLengthFromPoints=1;\n"
          <<"Mesh.ElementOrder=1;\n"
          <<"Mesh.SecondOrderIncomplete = 0;\n"
          <<"Mesh.Algorithm = 6;\n"
          <<"Mesh.OptimizeNetgen=1;\n"
          <<"// partitioning data\n"
          <<"Mesh.Partitioner=1;\n"
          <<"Mesh.NbPartitions=1;\n"
          <<"Mesh.MshFilePartitioned=0;\n"
          <<"h=" << h << ";\n"
          <<"Point(1) = {1,-1,0,h};\n"
          <<"Point(2) = {1,1,0,h};\n"
          <<"Point(3) = {-1,-1,0,h};\n"
          <<"Line(1) = {1,2};\n"
          <<"Line(2) = {2,3};\n"
          <<"Line(3) = {3,1};\n"
          <<"Transfinite Line{1} = 1;\n"
          <<"Transfinite Line{2} = 1;\n"
          <<"Transfinite Line{3} = 1;\n"
          <<"Line Loop(4) = {3,1,2};\n"
          <<"Plane Surface(5) = {4};\n"
          <<"Physical Line(\"vert\") = {1};\n"
          <<"Physical Line(\"hypo\") = {2};\n"
          <<"Physical Line(\"hor\") = {3};\n"
          <<"Physical Surface(9) = {5};\n";

    std::ostringstream nameStr;
    if ( std::abs( h - 2 ) < 1e-10 )
        nameStr << "one-elt-real-rot";
    else
        nameStr << "one-elt-mesh-rot";
    gmsh_ptrtype gmshp( new Gmsh );
    gmshp->setPrefix( nameStr.str() );
    gmshp->setDescription( costr.str() );
    return gmshp;
}


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
    po::options_description testhcurloptions("test h_curl options");
    testhcurloptions.add_options()
        ("hsize", po::value<double>()->default_value( 0.1 ), "mesh size")
        ("xmin", po::value<double>()->default_value( -1 ), "xmin of the reference element")
        ("ymin", po::value<double>()->default_value( -1 ), "ymin of the reference element")
        ("zmin", po::value<double>()->default_value( -1 ), "zmin of the reference element")
        ;
    return testhcurloptions.add( Feel::feel_options() );
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_hcurl" ,
                     "test_hcurl" ,
                     "0.1",
                     "Test for h_curl space (Dim=2 Order=1)",
                     AboutData::License_GPL,
                     "Copyright (c) 2009 Universite Joseph Fourier");
    about.addAuthor("Cecile Daversin", "developer", "cecile.daversin@lncmi.cnrs.fr", "");
    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;

}


class TestHCurl
    :
    public Application
{
    typedef Application super;

public:

    //! numerical type is double
    typedef double value_type;

    //! linear algebra backend factory
    typedef Backend<value_type> backend_type;
    //! linear algebra backend factory shared_ptr<> type
    typedef typename boost::shared_ptr<backend_type> backend_ptrtype ;

    //! geometry entities type composing the mesh, here Simplex in Dimension Dim of Order G_order
    typedef Simplex<2,1> convex_type;
    //! mesh type
    typedef Mesh<convex_type> mesh_type;
    //! mesh shared_ptr<> type
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    //! the basis type of our approximation space
    typedef bases<Nedelec<0> > basis_type;
    //! the approximation function space type
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    //! the approximation function space type (shared_ptr<> type)
    typedef boost::shared_ptr<space_type> space_ptrtype;
    //! an element type of the approximation function space
    typedef typename space_type::element_type element_type;

    //! the exporter factory type
    typedef Exporter<mesh_type> export_type;
    //! the exporter factory (shared_ptr<> type)
    typedef boost::shared_ptr<export_type> export_ptrtype;

    /**
     * Constructor
     */
    TestHCurl( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        M_backend( backend_type::build( this->vm() ) ),
        meshSize( this->vm()["hsize"].as<double>() ),
        exporter( Exporter<mesh_type>::New( this->vm() ) )
    {
        std::cout << "[TestHCurl]\n";

        this->changeRepository( boost::format( "%1%/h_%2%/" )
                                % this->about().appName()
                                % this->vm()["hsize"].as<double>()
                                );
    }

    /**
     * run the application
     */
    void shape_functions(gmsh_ptrtype (*one_element_mesh)(double));

private:
    //! linear algebra backend
    backend_ptrtype M_backend;

    //! mesh characteristic size
    double meshSize;

    //! exporter factory
    export_ptrtype exporter;

}; //TestHCurl

void
TestHCurl::shape_functions(gmsh_ptrtype (*one_element_mesh_desc_fun)(double))
{
    using namespace Feel::vf;

    mesh_ptrtype oneelement_mesh = createGMSHMesh( _mesh=new mesh_type,
                                                   _desc = one_element_mesh_desc_fun(2) );

    // then a fine mesh which we use to export the basis function to
    // visualize them
    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _desc=one_element_mesh_desc_fun(meshSize) );

    space_ptrtype Xh = space_type::New( oneelement_mesh ); // space associated with reference element

    std::cout << "Family = " << Xh->basis()->familyName() << "\n"
              << "Dim    = " << Xh->basis()->nDim << "\n"
              << "Order  = " << Xh->basis()->nOrder << "\n"
              << "NDof   = " << Xh->nLocalDof() << "\n";

    element_type U_ref( Xh, "U" );
    element_type V_ref(Xh, "V");

    // To store the shape functions
    // 0 : hypothenuse edge, 1 : vertical edge, 2 : horizontal edge
    std::vector<element_type> u_vec(3);

    std::string shape_name = "shape_functions";
    export_ptrtype exporter_shape( export_type::New( this->vm(),
                                                     (boost::format( "%1%-%2%-%3%" )
                                                      % this->about().appName()
                                                      % one_element_mesh_desc_fun(2)->prefix()
                                                      % shape_name).str() ) );

    exporter_shape->step(0)->setMesh( mesh );

    for( size_type i = 0;i < Xh->nLocalDof(); ++i )
        {
            // U_ref corresponds to shape function (on reference element)
            U_ref.zero();
            U_ref( i ) = 1;

            u_vec[i] = U_ref;

            std::ostringstream ostr;
            ostr <<  one_element_mesh_desc_fun(2)->prefix()<< "-" << Xh->basis()->familyName() << "-" << i;
            exporter_shape->step(0)->add( ostr.str(), U_ref );
        }

    exporter_shape->save();

    auto F = M_backend->newVector( Xh );

    //// *********************** Check  alpha_i(N_j) evaluations  on reference element (with idv keyword) ********////
    std::vector<double> checkidv(9);
    std::vector<double> checkform1(9);
    std::vector<std::string> edges = boost::assign::list_of("hypo")("vert")("hor");
    for( int i = 0;i < 3; ++i )
    {
        int edgeid = 0;
        BOOST_FOREACH( std::string edge, edges )
        {
            // on ref element
            auto v = integrate( markedfaces(oneelement_mesh, edge), trans(T())*(JinvT())*idv(u_vec[i])).evaluate()(0,0);
            if ( edgeid == i )
                BOOST_CHECK_CLOSE( v, 1, 1e-14 );
            else
                BOOST_CHECK_SMALL( v, 1e-14 );
            checkidv[3*i+edgeid] = v;
            form1( _test=Xh, _vector=F, _init=true) = integrate( markedfaces(oneelement_mesh, edge),
                                                                         trans(T())*(JinvT())*id(u_vec[i]));
            v = inner_product(u_vec[i], *F);
            if ( edgeid == i )
                BOOST_CHECK_CLOSE( v, 1, 1e-14 );
            else
                BOOST_CHECK_SMALL( v, 1e-14 );
            checkform1[3*i+edgeid] = v;

            ++edgeid;
        }

    }
    BOOST_TEST_MESSAGE( " ********** Values of alpha_i (N_j ) = delta_{i,j} (reference element) ********** \n"
                        << "\n"
                        << " ********** Using idv keyword ********************* "
                        << "\n");
        for( int i = 0;i < 3; ++i )
    {
        int edgeid = 0;
        BOOST_FOREACH( std::string edge, edges )
        {
            BOOST_TEST_MESSAGE( " *** dof N_"<< i << " (associated with " << edgeid << " edge) *** \n"
                                << "alpha_"<< edge << "(N_"<<i<<") = " << checkidv[3*i+edgeid] << "\n" );
            ++edgeid;
        }
        BOOST_TEST_MESSAGE( "*********************************************** \n" );
    }
    //// ************************************************************************************ ////
}


}
#if USE_BOOST_TEST

BOOST_AUTO_TEST_SUITE( space )

BOOST_AUTO_TEST_CASE( test_hcurl_N0_ref )
{
    BOOST_TEST_MESSAGE( "test_hcurl_N0 on reference element" );
    Feel::TestHCurl t(boost::unit_test::framework::master_test_suite().argc,
                      boost::unit_test::framework::master_test_suite().argv,
                      Feel::makeAbout(), Feel::makeOptions() );
    t.shape_functions(&Feel::oneelement_geometry_ref);
    BOOST_TEST_MESSAGE( "test_hcurl_N0 on reference element done" );
}
BOOST_AUTO_TEST_CASE( test_hcurl_N0_real )
{
    BOOST_TEST_MESSAGE( "test_hcurl_N0 on one real element" );
    Feel::TestHCurl t(boost::unit_test::framework::master_test_suite().argc,
                      boost::unit_test::framework::master_test_suite().argv,
                      Feel::makeAbout(), Feel::makeOptions() );
    t.shape_functions(&Feel::oneelement_geometry_real_1);
    BOOST_TEST_MESSAGE( "test_hcurl_N0 on one real element done" );
}

BOOST_AUTO_TEST_CASE( test_hcurl_N0_real_2 )
{
    BOOST_TEST_MESSAGE( "test_hcurl_N0 on one real element" );
    Feel::TestHCurl t(boost::unit_test::framework::master_test_suite().argc,
                      boost::unit_test::framework::master_test_suite().argv,
                      Feel::makeAbout(), Feel::makeOptions() );
    t.shape_functions(&Feel::oneelement_geometry_real_2);
    BOOST_TEST_MESSAGE( "test_hcurl_N0 on one real element done" );
}
BOOST_AUTO_TEST_SUITE_END()
#else

int
main( int argc, char* argv[] )
{

    Feel::TestHCurl app_hcurl( argc, argv, Feel::makeAbout(), Feel::makeOptions() );
    app_hcurl.shape_functions();
}

#endif
