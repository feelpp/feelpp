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
using namespace Feel;
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
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "test_hcurl" ,
                            "test_hcurl" ,
                            "0.1",
                           "Test for h_curl space (Dim=2 Order=1)",
                            Feel::AboutData::License_GPL,
                           "Copyright (c) 2009 Universite Joseph Fourier");
    about.addAuthor("Cecile Daversin", "developer", "cecile.daversin@lncmi.cnrs.fr", "");
    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;

}

class Test_Hcurl
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
    Test_Hcurl( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od ),
        M_backend( backend_type::build( this->vm() ) ),
        meshSize( this->vm()["hsize"].as<double>() ),
        exporter( Exporter<mesh_type>::New( this->vm() ) )
    {
        std::cout << "[Test_Hcurl]\n";

        this->changeRepository( boost::format( "%1%/h_%2%/" )
                                % this->about().appName()
                                % this->vm()["hsize"].as<double>()
                                );
    }

    /**
     * run the application
     */
    void shape_functions();
    gmsh_ptrtype oneelement_geometry_ref();
    gmsh_ptrtype oneelement_geometry_real_1();
    gmsh_ptrtype oneelement_geometry_real_2();

private:
    //! linear algebra backend
    backend_ptrtype M_backend;

    //! mesh characteristic size
    double meshSize;

    //! exporter factory
    export_ptrtype exporter;
 
}; //Test_Hcurl

void
Test_Hcurl::shape_functions()
{
    using namespace Feel::vf;
 
    mesh_ptrtype oneelement_mesh_ref = createGMSHMesh( _mesh=new mesh_type,
                                                       _desc = oneelement_geometry_ref());

    mesh_ptrtype oneelement_mesh_real_1 = createGMSHMesh( _mesh=new mesh_type,
                                                        _desc = oneelement_geometry_real_1());

    mesh_ptrtype oneelement_mesh_real_2 = createGMSHMesh( _mesh=new mesh_type,
                                                        _desc = oneelement_geometry_real_2());

    // then a fine mesh which we use to export the basis function to
    // visualize them
    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _desc=domain( _name="fine",
                                                      _shape="simplex",
                                                      _dim=2,
                                                      _h=meshSize,
                                                      _xmin=this->vm()["xmin"].as<double>(),
                                                      _ymin=this->vm()["ymin"].as<double>(),
                                                      _zmin=this->vm()["zmin"].as<double>() ) );

    space_ptrtype Xh_ref = space_type::New( oneelement_mesh_ref ); // space associated with reference element
    space_ptrtype Xh_real = space_type::New( oneelement_mesh_real_1 ); // space associated with real element 1 (homothetic transfo)
    space_ptrtype Xh_real2 = space_type::New( oneelement_mesh_real_2 ); // space associated with real element 2 (rotation)

    std::cout << "Family = " << Xh_ref->basis()->familyName() << "\n"
              << "Dim    = " << Xh_ref->basis()->nDim << "\n"
              << "Order  = " << Xh_ref->basis()->nOrder << "\n"
              << "NDof   = " << Xh_ref->nLocalDof() << "\n";

    element_type U_ref( Xh_ref, "U" );
    element_type V_ref(Xh_ref, "V");

    // To store the shape functions
    // 0 : hypothenuse edge, 1 : vertical edge, 2 : horizontal edge
    std::vector<element_type> u_vec(3);

    std::string shape_name = "shape_functions";
    export_ptrtype exporter_shape( export_type::New( this->vm(),
                                                     (boost::format( "%1%-%2%" )
                                                      % this->about().appName()
                                                      % shape_name).str() ) );
 
    exporter_shape->step(0)->setMesh( mesh );

    for( size_type i = 0;i < Xh_ref->nLocalDof(); ++i )
        {
            // U_ref corresponds to shape function (on reference element)
            U_ref.zero();
            U_ref( i ) = 1;

            u_vec[i] = U_ref;

            std::ostringstream ostr;
            ostr << Xh_ref->basis()->familyName() << "-" << i;
            exporter_shape->step(0)->add( ostr.str(), U_ref );
        }

    exporter_shape->save();



    std::cout << " ********** Test for shape functions (on reference element hat{K}) ********** \n"
              << "\n"
              << "********** Check for Jacobian and evaluation of shape function on middle of the edge **********"
              << "\n"
              << std::endl;

    //// *******************  Check each component of the integral (on reference element) (with idv keyword) ********************************* ////
    auto alpha0_N0_ref_print = 
        integrate( markedfaces(oneelement_mesh_ref, "hypo"), 
                   trans(vec(cst(-1.0/sqrt(2.0)), cst(1.0/sqrt(2.0))))*print(trans(J()), "jacobian transposed")
                   *print(idv(u_vec[0]), "N0 on hypo middle"), _Q<0>() ).evaluate();
    std::cout << "\n" << std::endl;

    auto alpha1_N0_ref_print = 
        integrate( markedfaces(oneelement_mesh_ref, "vert"), 
                   trans( vec(cst(0.0), cst(-1.0) ) )*print(trans(J()), "jacobian transposed")
                   *print(idv(u_vec[0]), "N0 on vert middle"), _Q<0>()).evaluate();
    std::cout << "\n" << std::endl;

    auto alpha2_N0_ref_print = 
        integrate( markedfaces(oneelement_mesh_ref, "hor"), 
                   trans( vec(cst(1.0), cst(0.0) ) )*print(trans(J()), "jacobian transposed")
                   *print(idv(u_vec[0]), "N0 on hor middle"), _Q<0>() ).evaluate();
    std::cout << "\n" << std::endl;

    auto alpha0_N1_ref_print = 
        integrate( markedfaces(oneelement_mesh_ref, "hypo"), 
                   trans(vec(cst(-1.0/sqrt(2.0)), cst(1.0/sqrt(2.0))))*print(trans(J()), "jacobian transposed")
                   *print(idv(u_vec[1]), "N1 on hypo middle"), _Q<0>() ).evaluate();
    std::cout << "\n" << std::endl;

    auto alpha1_N1_ref_print = 
        integrate( markedfaces(oneelement_mesh_ref, "vert"), 
                   trans( vec(cst(0.0), cst(-1.0) ) )*print(trans(J()), "jacobian transposed")
                   *print(idv(u_vec[1]), "N1 on vert middle"), _Q<0>()).evaluate();
    std::cout << "\n" << std::endl;

    auto alpha2_N1_ref_print = 
        integrate( markedfaces(oneelement_mesh_ref, "hor"), 
                   trans( vec(cst(1.0), cst(0.0) ) )*print(trans(J()), "jacobian transposed")
                   *print(idv(u_vec[1]), "N1 on hor middle"), _Q<0>() ).evaluate();
    std::cout << "\n" << std::endl;

    auto alpha0_N2_ref_print = 
        integrate( markedfaces(oneelement_mesh_ref, "hypo"), 
                   trans(vec(cst(-1.0/sqrt(2.0)), cst(1.0/sqrt(2.0))))*print(trans(J()), "jacobian transposed")
                   *print(idv(u_vec[2]), "N2 on hypo middle"), _Q<0>() ).evaluate();
    std::cout << "\n" << std::endl;

    auto alpha1_N2_ref_print = 
        integrate( markedfaces(oneelement_mesh_ref, "vert"), 
                   trans( vec(cst(0.0), cst(-1.0) ) )*print(trans(J()), "jacobian transposed")
                   *print(idv(u_vec[2]), "N2 on vert middle"), _Q<0>()).evaluate();
    std::cout << "\n" << std::endl;

    auto alpha2_N2_ref_print = 
        integrate( markedfaces(oneelement_mesh_ref, "hor"), 
                   trans( vec(cst(1.0), cst(0.0) ) )*print(trans(J()), "jacobian transposed")
                   *print(idv(u_vec[2]), "N2 on hor middle"), _Q<0>() ).evaluate();
    std::cout << "\n" << std::endl;
    //// ************************************************************************************ ////

    //// *********************** Check  alpha_i(N_j) evaluations  on reference element (with idv keyword) ********////
    auto alpha0_N0_ref_idv = 
        integrate( markedfaces(oneelement_mesh_ref, "hypo"), trans(vec(cst(-1.0/sqrt(2.0)), cst(1.0/sqrt(2.0))))*trans(J())*idv(u_vec[0])).evaluate();
    auto alpha1_N0_ref_idv = integrate( markedfaces(oneelement_mesh_ref, "vert"), trans( vec(cst(0.0), cst(-1.0) ) )*trans(J())*idv(u_vec[0])).evaluate();
    auto alpha2_N0_ref_idv = integrate( markedfaces(oneelement_mesh_ref, "hor"), trans( vec(cst(1.0), cst(0.0) ) )*trans(J())*idv(u_vec[0])).evaluate();

    auto alpha0_N1_ref_idv = 
        integrate( markedfaces(oneelement_mesh_ref, "hypo"), trans(vec(cst(-1.0/sqrt(2.0)), cst(1.0/sqrt(2.0))))*trans(J())*idv(u_vec[1])).evaluate();
    auto alpha1_N1_ref_idv = integrate( markedfaces(oneelement_mesh_ref, "vert"), trans( vec(cst(0.0), cst(-1.0) ) )*trans(J())*idv(u_vec[1])).evaluate();
    auto alpha2_N1_ref_idv = integrate( markedfaces(oneelement_mesh_ref, "hor"), trans( vec(cst(1.0), cst(0.0) ) )*trans(J())*idv(u_vec[1])).evaluate();

    auto alpha0_N2_ref_idv = 
        integrate( markedfaces(oneelement_mesh_ref, "hypo"), trans(vec(cst(-1.0/sqrt(2.0)), cst(1.0/sqrt(2.0))))*trans(J())*idv(u_vec[2])).evaluate();
    auto alpha1_N2_ref_idv = integrate( markedfaces(oneelement_mesh_ref, "vert"), trans( vec(cst(0.0), cst(-1.0) ) )*trans(J())*idv(u_vec[2])).evaluate();
    auto alpha2_N2_ref_idv = integrate( markedfaces(oneelement_mesh_ref, "hor"), trans( vec(cst(1.0), cst(0.0) ) )*trans(J())*idv(u_vec[2])).evaluate();

    std::cout << " ********** Values of alpha_i (N_j ) = delta_{i,j} (reference element) ********** \n"
              << "\n"
              << " ********** Using idv keyword ********************* "
              << "\n"
              << " *** dof N_0 (associated with hypotenuse edge) *** \n"
              << "alpha_0(N_0) = " << alpha0_N0_ref_idv << "\n"
              << "alpha_1(N_0) = " << alpha1_N0_ref_idv << "\n"
              << "alpha_2(N_0) = " << alpha2_N0_ref_idv << "\n"
              << "*********************************************** \n"
              << std::endl;

    std::cout << " *** dof N_1 (associated with vertical edge) *** \n"
              << "alpha_0(N_1) = " << alpha0_N1_ref_idv << "\n"
              << "alpha_1(N_1) = " << alpha1_N1_ref_idv << "\n"
              << "alpha_2(N_1) = " << alpha2_N1_ref_idv << "\n"
              << "*********************************************** \n"
              << std::endl;

    std::cout << " *** dof N_2 (associated with horizontal edge) *** \n"
              << "alpha_0(N_2) = " << alpha0_N2_ref_idv << "\n"
              << "alpha_1(N_2) = " << alpha1_N2_ref_idv << "\n"
              << "alpha_2(N_2) = " << alpha2_N2_ref_idv << "\n"
              << "*********************************************** \n"
              << std::endl;

    //// ************************************************************************************ ////

    //// *********************** Check each component on reference element (with form1 / id keywords) ********////

        std::cout << "Components of the integral in reference element (form1 / id keywords) \n " << std::endl;

    //shape function alpha0
    auto F_ref_0_print = M_backend->newVector( Xh_ref );
    form1( _test = Xh_ref, _vector = F_ref_0_print, _init=true) =
        integrate( markedfaces(oneelement_mesh_ref, "hypo"), trans(vec(cst(-1.0/sqrt(2.0)), cst(1.0/sqrt(2.0))))*trans(J())*
                   print(id(V_ref), "Ni on hypo middle"), _Q<0>() );
    F_ref_0_print->close();
    std::cout << "\n" << std::endl;

    //shape function alpha1
    auto F_ref_1_print = M_backend->newVector( Xh_ref );
    form1( _test = Xh_ref, _vector = F_ref_1_print, _init=true) =
        integrate( markedfaces(oneelement_mesh_ref, "vert"), trans(vec(cst(0.0), cst(-1.0)))*trans(J())*print(id(V_ref),"Ni on vert middle") , _Q<0>() );
    F_ref_1_print->close();
    std::cout << "\n" << std::endl;

    //shape function alpha2
    auto F_ref_2_print = M_backend->newVector( Xh_ref );
    form1( _test = Xh_ref, _vector = F_ref_2_print, _init=true) =
        integrate( markedfaces(oneelement_mesh_ref, "hor"), trans(vec(cst(1.0), cst(0.0)))*trans(J())*print(id(V_ref), "Ni on hor middle"), _Q<0>() );
    F_ref_2_print->close();
    std::cout << "\n" << std::endl;
    //// ************************************************************************************ ////

    //// *********************** Check  alpha_i(N_j) evaluations  on reference element (with form1 / id keywords) ********////
    //shape function alpha0
    auto F_ref_0 = M_backend->newVector( Xh_ref );
    form1( _test = Xh_ref, _vector = F_ref_0, _init=true) =
        integrate( markedfaces(oneelement_mesh_ref, "hypo"), trans(vec(cst(-1.0/sqrt(2.0)), cst(1.0/sqrt(2.0))))*trans(J())*id(V_ref) );
    F_ref_0->close();

    //shape function alpha1
    auto F_ref_1 = M_backend->newVector( Xh_ref );
    form1( _test = Xh_ref, _vector = F_ref_1, _init=true) =
        integrate( markedfaces(oneelement_mesh_ref, "vert"), trans(vec(cst(0.0), cst(-1.0)))*trans(J())*id(V_ref) );
    F_ref_1->close();

    //shape function alpha2
    auto F_ref_2 = M_backend->newVector( Xh_ref );
    form1( _test = Xh_ref, _vector = F_ref_2, _init=true) =
        integrate( markedfaces(oneelement_mesh_ref, "hor"), trans(vec(cst(1.0), cst(0.0)))*trans(J())*id(V_ref) );
    F_ref_2->close();

    auto alpha0_N0_ref_form1 = inner_product(u_vec[0], *F_ref_0);
    auto alpha1_N0_ref_form1 = inner_product(u_vec[0], *F_ref_1);
    auto alpha2_N0_ref_form1 = inner_product(u_vec[0], *F_ref_2);

    auto alpha0_N1_ref_form1 = inner_product(u_vec[1], *F_ref_0);
    auto alpha1_N1_ref_form1 = inner_product(u_vec[1], *F_ref_1);
    auto alpha2_N1_ref_form1 = inner_product(u_vec[1], *F_ref_2);

    auto alpha0_N2_ref_form1 = inner_product(u_vec[2], *F_ref_0);
    auto alpha1_N2_ref_form1 = inner_product(u_vec[2], *F_ref_1);
    auto alpha2_N2_ref_form1 = inner_product(u_vec[2], *F_ref_2);


    std::cout << " ********** Values of alpha_i (N_j ) = delta_{i,j} (reference element) ********** \n"
              << "\n"
              << " ********** Using form1 keyword ********************* "
              << "\n"
              << " *** dof N_0 (associated with hypotenuse edge) *** \n"
              << "alpha_0(N_0) = " << alpha0_N0_ref_form1 << "\n"
              << "alpha_1(N_0) = " << alpha1_N0_ref_form1 << "\n"
              << "alpha_2(N_0) = " << alpha2_N0_ref_form1 << "\n"
              << "*********************************************** \n"
              << std::endl;

    std::cout << " *** dof N_1 (associated with vertical edge) *** \n"
              << "alpha_0(N_1) = " << alpha0_N1_ref_form1 << "\n"
              << "alpha_1(N_1) = " << alpha1_N1_ref_form1 << "\n"
              << "alpha_2(N_1) = " << alpha2_N1_ref_form1 << "\n"
              << "*********************************************** \n"
              << std::endl;

    std::cout << " *** dof N_2 (associated with horizontal edge) *** \n"
              << "alpha_0(N_2) = " << alpha0_N2_ref_form1 << "\n"
              << "alpha_1(N_2) = " << alpha1_N2_ref_form1 << "\n"
              << "alpha_2(N_2) = " << alpha2_N2_ref_form1 << "\n"
              << "*********************************************** \n"
              << std::endl;

    std::cout << " ********** Test for shape functions (on real element 1 (homothetic transformation from hat{K} ) ********** \n"
              << "\n"
              << std::endl;

    //// ************************************************************************************ ////

    //// ********* Check components  in real element 1 (homothetic transformation) (idv keyword) *************** ////

    std::cout << "Components of the integral in real element 1 (homothetic transformation) (idv keyword) \n" << std::endl;

    auto alpha0_N0_real_print = 
        integrate( markedfaces(oneelement_mesh_real_1, "hypo"), 
                   trans(vec(cst(-1.0/sqrt(2.0)), cst(1.0/sqrt(2.0))))*print(trans(J()), "jacobian transposed")*
                   print( idv(u_vec[0]), "N0 middle hypo"), _Q<0>()).evaluate();
    std::cout << "\n" << std::endl;

    auto alpha1_N0_real_print = 
        integrate( markedfaces(oneelement_mesh_real_1, "vert"), 
                   trans( vec(cst(0.0), cst(-1.0) ) )*print(trans(J()), "jacobian transposed")*
                   print( idv(u_vec[0]), "N0 middle vert"), _Q<0>()).evaluate();
    std::cout << "\n" << std::endl;

    auto alpha2_N0_real_print =
        integrate( markedfaces(oneelement_mesh_real_1, "hor"),
                   trans( vec(cst(1.0), cst(0.0) ) )*print(trans(J()), "jacobian transposed")*
                   print(idv(u_vec[0]), "N0 middle hor"), _Q<0>()).evaluate();
    std::cout << "\n" << std::endl;

    auto alpha0_N1_real_print =
        integrate( markedfaces(oneelement_mesh_real_1, "hypo"),
                   trans(vec(cst(-1.0/sqrt(2.0)), cst(1.0/sqrt(2.0))))*print(trans(J()), "jacobian transposed")*
                   print(idv(u_vec[1]),"N1 middle hypo"), _Q<0>()).evaluate();
   std::cout << "\n" << std::endl;

    auto alpha1_N1_real_print =
        integrate( markedfaces(oneelement_mesh_real_1, "vert"),
                   trans( vec(cst(0.0), cst(-1.0) ) )*print(trans(J()), "jacobian transposed")*
                   print(idv(u_vec[1]),"N1 middle vert"), _Q<0>()).evaluate();
    std::cout << "\n" << std::endl;

    auto alpha2_N1_real_print =
        integrate( markedfaces(oneelement_mesh_real_1, "hor"),
                   trans( vec(cst(1.0), cst(0.0) ) )*print(trans(J()), "jacobian transposed")*
                   print(idv(u_vec[1]),"N1 middle hor"), _Q<0>()).evaluate();
    std::cout << "\n" << std::endl;

    auto alpha0_N2_real_print =
        integrate( markedfaces(oneelement_mesh_real_1, "hypo"),
                   trans(vec(cst(-1.0/sqrt(2.0)), cst(1.0/sqrt(2.0))))*print(trans(J()), "jacobian transposed")*
                   print(idv(u_vec[2]),"N2 middle hypo"), _Q<0>()).evaluate();
    std::cout << "\n" << std::endl;

    auto alpha1_N2_real_print =
        integrate( markedfaces(oneelement_mesh_real_1, "vert"),
                   trans( vec(cst(0.0), cst(-1.0) ) )*print(trans(J()), "jacobian transposed")*
                   print(idv(u_vec[2]),"N2 middle vert"), _Q<0>()).evaluate();
    std::cout << "\n" << std::endl;

    auto alpha2_N2_real_print =
        integrate( markedfaces(oneelement_mesh_real_1, "hor"),
                   trans( vec(cst(1.0), cst(0.0) ) )*print(trans(J()), "jacobian transposed")*
                   print(idv(u_vec[2]),"N2 middle hor"), _Q<0>()).evaluate();
    std::cout << "\n" << std::endl;

    //// ************************************************************************************ ////

    std::cout << "\n" << std::endl;

    //// *********** Check alpha_i(N_j) evaluations  on real element 1 (homothetic transformation) (with idv keyword) ************ ////
    auto alpha0_N0_real_idv = 
        integrate( markedfaces(oneelement_mesh_real_1, "hypo"), trans(vec(cst(-1.0/sqrt(2.0)), cst(1.0/sqrt(2.0))))*trans(J())*idv(u_vec[0])).evaluate();
    auto alpha1_N0_real_idv = integrate( markedfaces(oneelement_mesh_real_1, "vert"), trans( vec(cst(0.0), cst(-1.0) ) )*trans(J())*idv(u_vec[0])).evaluate();
    auto alpha2_N0_real_idv = integrate( markedfaces(oneelement_mesh_real_1, "hor"), trans( vec(cst(1.0), cst(0.0) ) )*trans(J())*idv(u_vec[0])).evaluate();

    auto alpha0_N1_real_idv = 
        integrate( markedfaces(oneelement_mesh_real_1, "hypo"), trans(vec(cst(-1.0/sqrt(2.0)), cst(1.0/sqrt(2.0))))*trans(J())*idv(u_vec[1])).evaluate();
    auto alpha1_N1_real_idv = integrate( markedfaces(oneelement_mesh_real_1, "vert"), trans( vec(cst(0.0), cst(-1.0) ) )*trans(J())*idv(u_vec[1])).evaluate();
    auto alpha2_N1_real_idv = integrate( markedfaces(oneelement_mesh_real_1, "hor"), trans( vec(cst(1.0), cst(0.0) ) )*trans(J())*idv(u_vec[1])).evaluate();

    auto alpha0_N2_real_idv = 
        integrate( markedfaces(oneelement_mesh_real_1, "hypo"), trans(vec(cst(-1.0/sqrt(2.0)), cst(1.0/sqrt(2.0))))*trans(J())*idv(u_vec[2])).evaluate();
    auto alpha1_N2_real_idv = integrate( markedfaces(oneelement_mesh_real_1, "vert"), trans( vec(cst(0.0), cst(-1.0) ) )*trans(J())*idv(u_vec[2])).evaluate();
    auto alpha2_N2_real_idv = integrate( markedfaces(oneelement_mesh_real_1, "hor"), trans( vec(cst(1.0), cst(0.0) ) )*trans(J())*idv(u_vec[2])).evaluate();

    std::cout << " ********** Values of alpha_i (N_j ) = delta_{i,j} (real element 1) ********** \n"
              << "\n"
              << " ********** Using idv keyword ********************* "
              << "\n"
              << " *** dof N_0 (associated with hypotenuse edge) *** \n"
              << "alpha_0(N_0) = " << alpha0_N0_real_idv << "\n"
              << "alpha_1(N_0) = " << alpha1_N0_real_idv << "\n"
              << "alpha_2(N_0) = " << alpha2_N0_real_idv << "\n"
              << "*********************************************** \n"
              << std::endl;

    std::cout << " *** dof N_1 (associated with vertical edge) *** \n"
              << "alpha_0(N_1) = " << alpha0_N1_real_idv << "\n"
              << "alpha_1(N_1) = " << alpha1_N1_real_idv << "\n"
              << "alpha_2(N_1) = " << alpha2_N1_real_idv << "\n"
              << "*********************************************** \n"
              << std::endl;

    std::cout << " *** dof N_2 (associated with horizontal edge) *** \n"
              << "alpha_0(N_2) = " << alpha0_N2_real_idv << "\n"
              << "alpha_1(N_2) = " << alpha1_N2_real_idv << "\n"
              << "alpha_2(N_2) = " << alpha2_N2_real_idv << "\n"
              << "*********************************************** \n"
              << std::endl;

    //// ************************************************************************************ ////

     //// ********* Check components  in real element 1 (homothetic transformation) (form1 / id keywords) *************** ////

    std::cout << "Components of the integral in real element 1 (homothetic transformation) (form1 / idv keyword) \n" << std::endl;

        //shape function alpha0
    auto F_real_0_print = M_backend->newVector( Xh_real );
    form1( _test = Xh_real, _vector = F_real_0_print, _init=true) =
        integrate( markedfaces(oneelement_mesh_real_1, "hypo"), trans(vec(cst(-1.0/sqrt(2.0)), cst(1.0/sqrt(2.0))))*
                   trans(J())*print(id(V_ref), "Ni middle hypo"), _Q<0>() );
    F_real_0_print->close();
    std::cout << "\n" << std::endl;

    //shape function alpha1
    auto F_real_1_print = M_backend->newVector( Xh_real );
    form1( _test = Xh_real, _vector = F_real_1_print, _init=true) =
        integrate( markedfaces(oneelement_mesh_real_1, "vert"), trans(vec(cst(0.0), cst(-1.0)))*trans(J())*print( id(V_ref), "Ni middle vert"), _Q<0>() );
    F_real_1_print->close();
   std::cout << "\n" << std::endl;

    //shape function alpha2
    auto F_real_2_print = M_backend->newVector( Xh_real);
    form1( _test = Xh_real, _vector = F_real_2_print, _init=true) =
        integrate( markedfaces(oneelement_mesh_real_1, "hor"), trans(vec(cst(1.0), cst(0.0)))*trans(J())*print( id(V_ref), "Ni middle hor"), _Q<0>() );
    F_real_2_print->close();
   std::cout << "\n" << std::endl;

    //// ************************************************************************************ ////

    //// *********** Check alpha_i(N_j) evaluations  on real element 1 (homothetic transformation) (with form1 / id keyword) ************ ////

    //shape function alpha0
    auto F_real_0 = M_backend->newVector( Xh_real );
    form1( _test = Xh_real, _vector = F_real_0, _init=true) =
        integrate( markedfaces(oneelement_mesh_real_1, "hypo"), trans(vec(cst(-1.0/sqrt(2.0)), cst(1.0/sqrt(2.0))))*trans(J())*id(V_ref) );
    F_real_0->close();

    //shape function alpha1
    auto F_real_1 = M_backend->newVector( Xh_real );
    form1( _test = Xh_real, _vector = F_real_1, _init=true) =
        integrate( markedfaces(oneelement_mesh_real_1, "vert"), trans(vec(cst(0.0), cst(-1.0)))*trans(J())*id(V_ref) );
    F_real_1->close();

    //shape function alpha2
    auto F_real_2 = M_backend->newVector( Xh_real );
    form1( _test = Xh_real, _vector = F_real_2, _init=true) =
        integrate( markedfaces(oneelement_mesh_real_1, "hor"), trans(vec(cst(1.0), cst(0.0)))*trans(J())*id(V_ref) );
    F_real_2->close();

    auto alpha0_N0_real_form1 = inner_product(u_vec[0], *F_real_0);
    auto alpha1_N0_real_form1 = inner_product(u_vec[0], *F_real_1);
    auto alpha2_N0_real_form1 = inner_product(u_vec[0], *F_real_2);

    auto alpha0_N1_real_form1 = inner_product(u_vec[1], *F_real_0);
    auto alpha1_N1_real_form1 = inner_product(u_vec[1], *F_real_1);
    auto alpha2_N1_real_form1 = inner_product(u_vec[1], *F_real_2);

    auto alpha0_N2_real_form1 = inner_product(u_vec[2], *F_real_0);
    auto alpha1_N2_real_form1 = inner_product(u_vec[2], *F_real_1);
    auto alpha2_N2_real_form1 = inner_product(u_vec[2], *F_real_2);

    std::cout << " ********** Values of alpha_i (N_j ) = delta_{i,j} (real element 1) ********** \n"
              << "\n"
              << " ********** Using form1 keyword ********************* "
              << "\n"
              << " *** dof N_0 (associated with hypotenuse edge) *** \n"
              << "alpha_0(N_0) = " << alpha0_N0_real_form1 << "\n"
              << "alpha_1(N_0) = " << alpha1_N0_real_form1 << "\n"
              << "alpha_2(N_0) = " << alpha2_N0_real_form1 << "\n"
              << "*********************************************** \n"
              << std::endl;

    std::cout << " *** dof N_1 (associated with vertical edge) *** \n"
              << "alpha_0(N_1) = " << alpha0_N1_real_form1 << "\n"
              << "alpha_1(N_1) = " << alpha1_N1_real_form1 << "\n"
              << "alpha_2(N_1) = " << alpha2_N1_real_form1 << "\n"
              << "*********************************************** \n"
              << std::endl;

    std::cout << " *** dof N_2 (associated with horizontal edge) *** \n"
              << "alpha_0(N_2) = " << alpha0_N2_real_form1 << "\n"
              << "alpha_1(N_2) = " << alpha1_N2_real_form1 << "\n"
              << "alpha_2(N_2) = " << alpha2_N2_real_form1 << "\n"
              << "*********************************************** \n"
              << std::endl;

    //// ************************************************************************************ ////

    std::cout << " ********** Test for shape functions (on real element 2 (rotation pi/2 from hat{K} ) ********** \n"
              << "\n"
              << std::endl;

    //// ************************************************************************************ ////

    //// ********* Check components  in real element 2 (rotation pi/2) (idv keyword) *************** ////

    std::cout << "Components of the integral in real element 2 (rotation pi/2) (idv keyword) \n" << std::endl;

    auto alpha0_N0_real2_print =
        integrate( markedfaces(oneelement_mesh_real_2, "hypo"),
                   trans(vec(cst(-1.0/sqrt(2.0)), cst(-1.0/sqrt(2.0))))*print(trans(J()), "jacobian transposed")*
                   print( idv(u_vec[0]), "N0 middle hypo"), _Q<0>()).evaluate();
    std::cout << "\n" << std::endl;

    auto alpha1_N0_real2_print = 
        integrate( markedfaces(oneelement_mesh_real_2, "vert"), 
                   trans( vec(cst(0.0), cst(1.0) ) )*print(trans(J()), "jacobian transposed")*
                   print( idv(u_vec[0]), "N0 middle vert"), _Q<0>()).evaluate();
    std::cout << "\n" << std::endl;

    auto alpha2_N0_real2_print =
        integrate( markedfaces(oneelement_mesh_real_2, "hor"),
                   trans( vec(cst(1.0), cst(0.0) ) )*print(trans(J()), "jacobian transposed")*
                   print(idv(u_vec[0]), "N0 middle hor"), _Q<0>()).evaluate();
    std::cout << "\n" << std::endl;

    auto alpha0_N1_real2_print =
        integrate( markedfaces(oneelement_mesh_real_2, "hypo"),
                   trans(vec(cst(-1.0/sqrt(2.0)), cst(-1.0/sqrt(2.0))))*print(trans(J()), "jacobian transposed")*
                   print(idv(u_vec[1]),"N1 middle hypo"), _Q<0>()).evaluate();
   std::cout << "\n" << std::endl;

    auto alpha1_N1_real2_print =
        integrate( markedfaces(oneelement_mesh_real_2, "vert"),
                   trans( vec(cst(0.0), cst(1.0) ) )*print(trans(J()), "jacobian transposed")*
                   print(idv(u_vec[1]),"N1 middle vert"), _Q<0>()).evaluate();
    std::cout << "\n" << std::endl;

    auto alpha2_N1_real2_print =
        integrate( markedfaces(oneelement_mesh_real_2, "hor"),
                   trans( vec(cst(1.0), cst(0.0) ) )*print(trans(J()), "jacobian transposed")*
                   print(idv(u_vec[1]),"N1 middle hor"), _Q<0>()).evaluate();
    std::cout << "\n" << std::endl;

    auto alpha0_N2_real2_print =
        integrate( markedfaces(oneelement_mesh_real_2, "hypo"),
                   trans(vec(cst(-1.0/sqrt(2.0)), cst(-1.0/sqrt(2.0))))*print(trans(J()), "jacobian transposed")*
                   print(idv(u_vec[2]),"N2 middle hypo"), _Q<0>()).evaluate();
    std::cout << "\n" << std::endl;

    auto alpha1_N2_real2_print =
        integrate( markedfaces(oneelement_mesh_real_2, "vert"),
                   trans( vec(cst(0.0), cst(1.0) ) )*print(trans(J()), "jacobian transposed")*
                   print(idv(u_vec[2]),"N2 middle vert"), _Q<0>()).evaluate();
    std::cout << "\n" << std::endl;

    auto alpha2_N2_real2_print =
        integrate( markedfaces(oneelement_mesh_real_2, "hor"),
                   trans( vec(cst(1.0), cst(0.0) ) )*print(trans(J()), "jacobian transposed")*
                   print(idv(u_vec[2]),"N2 middle hor"), _Q<0>()).evaluate();
    std::cout << "\n" << std::endl;

    //// ************************************************************************************ ////

    std::cout << "\n" << std::endl;

    //// *********** Check alpha_i(N_j) evaluations  on real element 2 (rotation pi/2) (with idv keyword) ************ ////
    auto alpha0_N0_real2_idv = 
        integrate( markedfaces(oneelement_mesh_real_2, "hypo"), trans(vec(cst(-1.0/sqrt(2.0)), cst(-1.0/sqrt(2.0))))*trans(J())*idv(u_vec[0])).evaluate();
    auto alpha1_N0_real2_idv = integrate( markedfaces(oneelement_mesh_real_2, "vert"), trans( vec(cst(0.0), cst(1.0) ) )*trans(J())*idv(u_vec[0])).evaluate();
    auto alpha2_N0_real2_idv = integrate( markedfaces(oneelement_mesh_real_2, "hor"), trans( vec(cst(1.0), cst(0.0) ) )*trans(J())*idv(u_vec[0])).evaluate();

    auto alpha0_N1_real2_idv = 
        integrate( markedfaces(oneelement_mesh_real_2, "hypo"), trans(vec(cst(-1.0/sqrt(2.0)), cst(-1.0/sqrt(2.0))))*trans(J())*idv(u_vec[1])).evaluate();
    auto alpha1_N1_real2_idv = integrate( markedfaces(oneelement_mesh_real_2, "vert"), trans( vec(cst(0.0), cst(1.0) ) )*trans(J())*idv(u_vec[1])).evaluate();
    auto alpha2_N1_real2_idv = integrate( markedfaces(oneelement_mesh_real_2, "hor"), trans( vec(cst(1.0), cst(0.0) ) )*trans(J())*idv(u_vec[1])).evaluate();

    auto alpha0_N2_real2_idv = 
        integrate( markedfaces(oneelement_mesh_real_2, "hypo"), trans(vec(cst(-1.0/sqrt(2.0)), cst(-1.0/sqrt(2.0))))*trans(J())*idv(u_vec[2])).evaluate();
    auto alpha1_N2_real2_idv = integrate( markedfaces(oneelement_mesh_real_2, "vert"), trans( vec(cst(0.0), cst(1.0) ) )*trans(J())*idv(u_vec[2])).evaluate();
    auto alpha2_N2_real2_idv = integrate( markedfaces(oneelement_mesh_real_2, "hor"), trans( vec(cst(1.0), cst(0.0) ) )*trans(J())*idv(u_vec[2])).evaluate();

    std::cout << " ********** Values of alpha_i (N_j ) = delta_{i,j} (real element 2) ********** \n"
              << "\n"
              << " ********** Using idv keyword ********************* "
              << "\n"
              << " *** dof N_0 (associated with hypotenuse edge) *** \n"
              << "alpha_0(N_0) = " << alpha0_N0_real2_idv << "\n"
              << "alpha_1(N_0) = " << alpha1_N0_real2_idv << "\n"
              << "alpha_2(N_0) = " << alpha2_N0_real2_idv << "\n"
              << "*********************************************** \n"
              << std::endl;

    std::cout << " *** dof N_1 (associated with vertical edge) *** \n"
              << "alpha_0(N_1) = " << alpha0_N1_real2_idv << "\n"
              << "alpha_1(N_1) = " << alpha1_N1_real2_idv << "\n"
              << "alpha_2(N_1) = " << alpha2_N1_real2_idv << "\n"
              << "*********************************************** \n"
              << std::endl;

    std::cout << " *** dof N_2 (associated with horizontal edge) *** \n"
              << "alpha_0(N_2) = " << alpha0_N2_real2_idv << "\n"
              << "alpha_1(N_2) = " << alpha1_N2_real2_idv << "\n"
              << "alpha_2(N_2) = " << alpha2_N2_real2_idv << "\n"
              << "*********************************************** \n"
              << std::endl;

    //// ************************************************************************************ ////

     //// ********* Check components  in real element 2 (rotation pi/2) (form1 / id keywords) *************** ////

    std::cout << "Components of the integral in real element 1 (rotation pi/2) (form1 / idv keyword) \n" << std::endl;

        //shape function alpha0
    auto F_real2_0_print = M_backend->newVector( Xh_real2 );
    form1( _test = Xh_real2, _vector = F_real2_0_print, _init=true) =
        integrate( markedfaces(oneelement_mesh_real_2, "hypo"), trans(vec(cst(-1.0/sqrt(2.0)), cst(-1.0/sqrt(2.0))))*
                   trans(J())*print(id(V_ref), "Ni middle hypo"), _Q<0>() );
    F_real2_0_print->close();
    std::cout << "\n" << std::endl;

    //shape function alpha1
    auto F_real2_1_print = M_backend->newVector( Xh_real2 );
    form1( _test = Xh_real2, _vector = F_real2_1_print, _init=true) =
        integrate( markedfaces(oneelement_mesh_real_2, "vert"), trans(vec(cst(0.0), cst(1.0)))*trans(J())*print( id(V_ref), "Ni middle vert"), _Q<0>() );
    F_real2_1_print->close();
   std::cout << "\n" << std::endl;

    //shape function alpha2
    auto F_real2_2_print = M_backend->newVector( Xh_real2 );
    form1( _test = Xh_real2, _vector = F_real2_2_print, _init=true) =
        integrate( markedfaces(oneelement_mesh_real_2, "hor"), trans(vec(cst(1.0), cst(0.0)))*trans(J())*print( id(V_ref), "Ni middle hor"), _Q<0>() );
    F_real2_2_print->close();
   std::cout << "\n" << std::endl;

    //// ************************************************************************************ ////

    //// *********** Check alpha_i(N_j) evaluations  on real element 2 (rotation pi/2) (with form1 / id keyword) ************ ////

    //shape function alpha0
    auto F_real2_0 = M_backend->newVector( Xh_real2 );
    form1( _test = Xh_real2, _vector = F_real2_0, _init=true) =
        integrate( markedfaces(oneelement_mesh_real_2, "hypo"), trans(vec(cst(-1.0/sqrt(2.0)), cst(-1.0/sqrt(2.0))))*trans(J())*id(V_ref) );
    F_real2_0->close();

    //shape function alpha1
    auto F_real2_1 = M_backend->newVector( Xh_real2 );
    form1( _test = Xh_real2, _vector = F_real2_1, _init=true) =
        integrate( markedfaces(oneelement_mesh_real_2, "vert"), trans(vec(cst(0.0), cst(1.0)))*trans(J())*id(V_ref) );
    F_real2_1->close();

    //shape function alpha2
    auto F_real2_2 = M_backend->newVector( Xh_real2 );
    form1( _test = Xh_real2, _vector = F_real2_2, _init=true) =
        integrate( markedfaces(oneelement_mesh_real_2, "hor"), trans(vec(cst(1.0), cst(0.0)))*trans(J())*id(V_ref) );
    F_real2_2->close();

    auto alpha0_N0_real2_form1 = inner_product(u_vec[0], *F_real2_0);
    auto alpha1_N0_real2_form1 = inner_product(u_vec[0], *F_real2_1);
    auto alpha2_N0_real2_form1 = inner_product(u_vec[0], *F_real2_2);

    auto alpha0_N1_real2_form1 = inner_product(u_vec[1], *F_real2_0);
    auto alpha1_N1_real2_form1 = inner_product(u_vec[1], *F_real2_1);
    auto alpha2_N1_real2_form1 = inner_product(u_vec[1], *F_real2_2);

    auto alpha0_N2_real2_form1 = inner_product(u_vec[2], *F_real2_0);
    auto alpha1_N2_real2_form1 = inner_product(u_vec[2], *F_real2_1);
    auto alpha2_N2_real2_form1 = inner_product(u_vec[2], *F_real2_2);

    std::cout << " ********** Values of alpha_i (N_j ) = delta_{i,j} (real element 2) ********** \n"
              << "\n"
              << " ********** Using form1 keyword ********************* "
              << "\n"
              << " *** dof N_0 (associated with hypotenuse edge) *** \n"
              << "alpha_0(N_0) = " << alpha0_N0_real2_form1 << "\n"
              << "alpha_1(N_0) = " << alpha1_N0_real2_form1 << "\n"
              << "alpha_2(N_0) = " << alpha2_N0_real2_form1 << "\n"
              << "*********************************************** \n"
              << std::endl;

    std::cout << " *** dof N_1 (associated with vertical edge) *** \n"
              << "alpha_0(N_1) = " << alpha0_N1_real2_form1 << "\n"
              << "alpha_1(N_1) = " << alpha1_N1_real2_form1 << "\n"
              << "alpha_2(N_1) = " << alpha2_N1_real2_form1 << "\n"
              << "*********************************************** \n"
              << std::endl;

    std::cout << " *** dof N_2 (associated with horizontal edge) *** \n"
              << "alpha_0(N_2) = " << alpha0_N2_real2_form1 << "\n"
              << "alpha_1(N_2) = " << alpha1_N2_real2_form1 << "\n"
              << "alpha_2(N_2) = " << alpha2_N2_real2_form1 << "\n"
              << "*********************************************** \n"
              << std::endl;

    //// ************************************************************************************ ////
}


/// Geometry for one-element meshes
gmsh_ptrtype
Test_Hcurl::oneelement_geometry_ref()
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
          <<"h=2;\n"
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
    nameStr << "one-elt-ref";
    gmsh_ptrtype gmshp( new Gmsh );
    gmshp->setPrefix( nameStr.str() );
    gmshp->setDescription( costr.str() );
    return gmshp;
}

// homothetic transformation of reference element (center 0, rate 2)
gmsh_ptrtype
Test_Hcurl::oneelement_geometry_real_1()
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
          <<"h=2;\n"
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
    nameStr << "one-elt-real-homo";
    gmsh_ptrtype gmshp( new Gmsh );
    gmshp->setPrefix( nameStr.str() );
    gmshp->setDescription( costr.str() );
    return gmshp;
}

// Rotation of angle (pi/2)
gmsh_ptrtype
Test_Hcurl::oneelement_geometry_real_2()
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
          <<"h=2;\n"
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
    nameStr << "one-elt-real-rot";
    gmsh_ptrtype gmshp( new Gmsh );
    gmshp->setPrefix( nameStr.str() );
    gmshp->setDescription( costr.str() );
    return gmshp;
}

int
main( int argc, char* argv[] )
{

    Test_Hcurl app_hcurl( argc, argv, makeAbout(), makeOptions() );
    app_hcurl.shape_functions();
}
