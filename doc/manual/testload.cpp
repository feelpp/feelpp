/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   Date: 2011-06-21

   Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
 \file testload.cpp
 \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
 \date 16-06-2011
 */
#include <feel/feelcore/feel.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelvf/vf.hpp>

int main(int argc, char** argv)
{
    double hsize = 2;
    // Declare the supported options.
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("filename", po::value<std::string>()->default_value( "Cylref.mesh" ), "h size")
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    using namespace Feel;
    using namespace Feel::vf;
    Feel::Environment env(argc, argv );
    typedef Mesh<Simplex<3> > mesh_type;
    std::string mesh_name=vm["filename"].as<std::string>();
    auto mesh = createGMSHMesh( _mesh=new mesh_type,
                                _desc=mshconvert( _filename=mesh_name ),
                                _physical_are_elementary_regions=true,
                                _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES );
    std::cout << "mesh " << mesh_name << " loaded\n" << std::endl;

    std::cout << "volume =" << std::endl
              << integrate( elements(mesh), cst(1.) ).evaluate() << "\n";
    std::cout << "surface =" << std::endl
              << integrate( boundaryfaces(mesh), cst(1.) ).evaluate() << "\n";


}
#if 0
#include <feel/options.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/application.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feelfilters/gmsh.hpp>

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
    po::options_description myappoptions("testload options");
    return myappoptions.add( feel_options() ).add( backend_options( "testload" ) );
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
    AboutData about( "testload" ,
					"testload" ,
					"0.1",
					"test to load medit meshes",
					AboutData::License_GPL,
					"Copyright (c) 2011 Universite Joseph Fourier");

    about.addAuthor("Christophe Prud'homme",
                    "developer",
                    "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;
}

/**
 * \class testload
 *
 * Test for loading a medit mesh
 */
template<int Dim>
class TestLoad: public Application //Feel::Simget
{
public:

	typedef Simplex<Dim> convex_type;
	typedef Mesh<convex_type > mesh_type;
	typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    /**
     * constructor only about data and no options description
     */
    TestLoad( int argc, char** argv, AboutData const& );

    /**
     * constructor about data and options description
     */
    TestLoad( int argc, char** argv, AboutData const&,
						po::options_description const&  );

    /**
     * This function is responsible for the actual work done by TestLoad.
     */
    void run();

};
template<int Dim>
TestLoad<Dim>::TestLoad(int argc, char** argv,
             AboutData const& ad )
:
Application( argc, argv, ad )
{}

template<int Dim>
TestLoad<Dim>::TestLoad(int argc, char** argv,
             AboutData const& ad,
             po::options_description const& od )
:
Application( argc, argv, ad, od )
{}

template<int Dim>
void TestLoad<Dim>::run()
{
    /**
     * print the help if --help is passed as an argument
     */
    if ( this->vm().count( "help" ) )
	{
		std::cout << this->optionsDescription() << "\n";
		return;
	}

    /**
     * store all subsequent data files in a HOME/feel/doc/tutorial/TestLoad/
     */
    this->changeRepository( boost::format( "doc/tutorial/%1%/" )
						   % this->about().appName() );


	/**
	 * The mesh to be loaded must be placed in ~/feel/doc/tutorial/testload/
	 */
	mesh_ptrtype mesh = loadGMSHMesh( _mesh=new mesh_type,
                                      _filename="Cylref.mesh",
                                      _physical_are_elementary_regions=true);

	/**
	 * Test for loading a classical GMSH since the changes
	 */
    //mesh_ptrtype mesh = loadGMSHMesh( _mesh= new mesh_type,
    //                                  _filename="hypercube-2.msh");

	Log() << "MESH LOADED\n";
}


/**
 * main function: entry point of the program
 */
int main( int argc, char** argv )
{
    Feel::Environment env( argc, argv );

    /**
     * intantiate a TestLoad class
     */
     TestLoad<3> app( argc, argv, makeAbout());

    /**
     * run the application
     */
    app.run();
}
#endif
