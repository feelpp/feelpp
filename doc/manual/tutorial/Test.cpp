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
   \file Test.cpp
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
    Feel::po::options_description TestOptions( "Test options" );
    TestOptions.add_options()
    ( "hsize", Feel::po::value<double>()->default_value( 0.1 ), "first h value to start convergence" )
    ( "GOrder", Feel::po::value<double>()->default_value( 1 ), "first geometrical Order value to start convergence" )
    ( "export-matlab", "export matrix and vectors in matlab" )
    ;
    return TestOptions.add( Feel::feel_options() ) ;
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
    Feel::AboutData about( "Test" ,
                           "Test" ,
                           "0.1",
                           "Test on the convergence of the geometrical order",
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
class Test
    :
public Application
{
    typedef Application super;
public:

    typedef double value_type;
    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;


    /*mesh*/
    typedef Simplex<2,4> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    
    /* export */
    typedef Exporter<mesh_type> export_type;

    FEELPP_DONT_INLINE
    Test( int argc, char** argv, AboutData const& ad, po::options_description const& od );

    // init mesh and space
    FEELPP_DONT_INLINE
    void init();

    /**
     * run the convergence test
     */
    FEELPP_DONT_INLINE
    void run();

private:

    double meshSize;
    mesh_ptrtype mesh;

    boost::shared_ptr<export_type> exporter;
}; // Stokes


Test::Test( int argc, char** argv, AboutData const& ad, po::options_description const& od )
    :
    super( argc, argv, ad, od ),
    meshSize( this->vm()["hsize"].as<double>() ),
    exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) )
{
}

void
Test::init()
{
    if ( this->vm().count( "help" ) )
    {
        std::cout << this->optionsDescription() << "\n";
        return;
    }

    if ( this->vm().count( "nochdir" ) == false )
        this->changeRepository( boost::format( "doc/tutorial/%1%/%2%/Part%4%/h_%3%/" )
                                % this->about().appName()
                                % convex_type::name()
                                % this->vm()["hsize"].as<double>()
                                % Environment::numberOfProcessors() );

    //************************ Cylindre *************************************** 
    // GeoTool::Node Centre(0,0,0);
    // GeoTool::Node Rayon( 1);
    // GeoTool::Node Dir(1,0,0);
    // GeoTool::Node Lg(5,0,0);
    // GeoTool::Cylindre C( meshSize,"Cyl",Centre,Dir,Rayon,Lg);
    // C.setMarker(_type="surface",_name="Inlet",_marker1=true);
    // C.setMarker(_type="surface",_name="Outlet",_marker2=true);
    // C.setMarker(_type="surface",_name="Wall",_marker3=true);
    // C.setMarker(_type="volume",_name="OmegaFluide",_markerAll=true);

    // mesh = C.createMesh(_mesh= new mesh_type,
    //                     _name="cylinder",
    //                     _partitions=Environment::worldComm().localSize(),
    //                     _worldcomm=Environment::worldComm() );
    //*************************************************************************


    //************************ Cercle ***************************************
    GeoTool::Node xc(0,0);
    GeoTool::Node xr(3,0);

    GeoTool::Circle Circ(meshSize,"OmegaCircle",xc,xr);

    Circ.setMarker(_type="line",_name="Boundary2",_markerAll=true);
    Circ.setMarker(_type="surface",_name="Omega",_markerAll=true);

    //xc est le centre et xr est un point sur le cercle
    mesh = Circ.createMesh(_mesh= new mesh_type,
                           _name="circle",
                           _partitions=Environment::worldComm().localSize(),
                           _worldcomm=Environment::worldComm() );
    //*************************************************************************


}
void
Test::run()
{
    this->init();

    //cylinder
    // double exact_volume = 5*4*math::atan(1.);    //pi
    // std::cout << "Vexact=" << exact_volume<< "\n";
    // double local_volume =  integrate( elements( mesh ), constant( 1.0 ) ).evaluate( false )( 0,0 );
    // double global_volume = integrate( elements( mesh ), constant( 1.0 ) ).evaluate()( 0,0 );
    // std::cout << "Global volume= " << global_volume<< "[ " << local_volume << " ]\n";

    // auto diff = exact_volume - global_volume;

    // std::cout << "||V-Vexact||=" <<diff<< "\n";


    //disque
    double exact_surface = 9*4*math::atan(1.);    //pi
    std::cout <<"\n"<<"Sexact=" << exact_surface<< "\n";
    double local_surface =  integrate( elements( mesh ), constant( 1.0 ) ).evaluate( false )( 0,0 );
    double global_surface = integrate( elements( mesh ), constant( 1.0 ) ).evaluate()( 0,0 );
    // double global_surface = integrate( elements( mesh ), constant( 1.0 ),_quad=_Q<20>()  ).evaluate()( 0,0 );
    std::cout << "Global surface= " << global_surface<< "[ " << local_surface << " ]\n";

    auto diff = exact_surface - global_surface;

    std::cout << "||S-Sexact||=" <<diff<< "\n";

} // Stokes::run
}

int
main( int argc, char** argv )
{

    using namespace Feel;
    Environment env( argc, argv );

    /* assertions handling */
    Feel::Assert::setLog( "Test.assert" );

    /* define and run application */
    Feel::Test Test( argc, argv, makeAbout(), makeOptions() );
    Test.run();
}
