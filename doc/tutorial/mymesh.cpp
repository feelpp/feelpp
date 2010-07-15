/* -*- mode: c++ coding: utf-8 -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2008-02-06

  Copyright (C) 2008-2009 Université Joseph Fourier (Grenoble I)


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
   \file mymesh.cpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2008-02-06
 */
#include <life/options.hpp>
#include <life/lifecore/application.hpp>
#include <life/lifediscr/mesh.hpp>
#include <life/lifediscr/functionspace.hpp>
#include <life/lifediscr/region.hpp>
#include <life/lifefilters/gmsh.hpp>

#include <life/lifefilters/exporter.hpp>

using namespace Life;

inline
Life::po::options_description
makeOptions()
{
    Life::po::options_description mymeshoptions("MyMesh options");
    mymeshoptions.add_options()
        ("hsize", Life::po::value<double>()->default_value( 0.1 ), "mesh size")
        ("shape", Life::po::value<std::string>()->default_value( "hypercube" ), "shape of the domain (either simplex or hypercube)")
        ;

    // return the options mymeshoptions and the life_options defined
    // internally by Life
    return mymeshoptions.add( Life::life_options() );
}
inline
Life::AboutData
makeAbout()
{
    Life::AboutData about( "mymesh" ,
                           "mymesh" ,
                           "0.2",
                           "my first Life application",
                           Life::AboutData::License_GPL,
                           "Copyright (c) 2008-2010 Universite Joseph Fourier");

    about.addAuthor("Christophe Prud'homme",
                    "developer",
                    "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;
}

template<int Dim>
class MyMesh: public Life::Simget
{
public:

    //# marker1 #
    typedef Simplex<Dim> convex_type;
    //typedef SimplexProduct<Dim, 1,Dim> convex_type;
    //# endmarker1 #

    //# marker2 #
    typedef Mesh<convex_type > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    //# endmarker2 #

    //# marker61 #
    /* export */
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;
    //# endmarker61 #

    /**
     * constructor about data and options description
     */
    MyMesh( po::variables_map const& vm, AboutData const& about  );

    void run();

    void run( const double* X, unsigned long P, double* Y, unsigned long N );

private:
    double meshSize;
    std::string shape;
    export_ptrtype exporter;
};

template<int Dim>
MyMesh<Dim>::MyMesh( po::variables_map const& vm, AboutData const& about )
    :
    Simget( vm, about ),
    meshSize( this->vm()["hsize"].template as<double>() ),
    shape( this->vm()["shape"].template as<std::string>() ),
    exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) )
{}


template<int Dim>
void
MyMesh<Dim>::run()
{
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Execute MyMesh<" << Dim << ">\n";
    std::vector<double> X( 2 );
    X[0] = meshSize;
    if ( shape == "hypercube" )
        X[1] = 1;
    else // default is simplex
        X[1] = 0;
    std::vector<double> Y( 3 );
    run( X.data(), X.size(), Y.data(), Y.size() );
}
template<int Dim>
void
MyMesh<Dim>::run( const double* X, unsigned long P, double* Y, unsigned long N )
{
    if ( X[1] == 0 ) shape = "simplex";
    if ( X[1] == 1 ) shape = "hypercube";

    if ( !this->vm().count( "nochdir" ) )
        Environment::changeRepository( boost::format( "doc/tutorial/%1%/%2%-%3%/h_%4%/" )
                                       % this->about().appName()
                                       % shape
                                       % Dim
                                       % meshSize );
    //# marker4 #
    mesh_ptrtype mesh = createGMSHMesh( _mesh=new mesh_type,
                                        _desc=domain( _name=(boost::format( "%1%-%2%" ) % shape % Dim).str() ,
                                                      _shape=shape,
                                                      _dim=Dim,
                                                      _h=X[0] ) );
    mesh->setComponents( MESH_PARTITION| MESH_UPDATE_FACES|MESH_UPDATE_EDGES);
    mesh->updateForUse();
    //#endmarker4#


    //# marker62 #
    exporter->step(0)->setMesh( mesh );
    exporter->save();
    //# endmarker62 #
}

//
// main function: entry point of the program
//
int main( int argc, char** argv )
{
    Application app( argc, argv, makeAbout(), makeOptions() );
    if ( app.vm().count( "help" ) )
    {
        std::cout << app.optionsDescription() << "\n";
        return 0;
    }

    app.add( new MyMesh<1>( app.vm(), app.about() ) );
    app.add( new MyMesh<2>( app.vm(), app.about() ) );
    app.add( new MyMesh<3>( app.vm(), app.about() ) );

    app.run();
}


