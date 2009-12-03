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
#include <life/lifefilters/gmshtensorizeddomain.hpp>
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
                           "0.1",
                           "my first Life application",
                           Life::AboutData::License_GPL,
                           "Copyright (c) 2008 Universite Joseph Fourier");

    about.addAuthor("Christophe Prud'homme",
                    "developer",
                    "christophe.prudhomme@ujf-grenoble.fr", "");
    return about;
}

template<int Dim>
class MyMesh: public Life::Application
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

    /* export */
    typedef Exporter<mesh_type> export_type;
    typedef boost::shared_ptr<export_type> export_ptrtype;

    /**
     * constructor about data and options description
     */
    MyMesh( int argc, char** argv,
           Life::AboutData const&,
           Life::po::options_description const&  );

    /**
     * create the mesh
     */
    mesh_ptrtype createMesh();

    /**
     * must be redefined by Application subclass
     */
    void run();

private:

    /**
     * set the application member data from command line options
     */
    void setOptions();

private:
    double M_meshSize;

    export_ptrtype exporter;
};

template<int Dim>
MyMesh<Dim>::MyMesh(int argc, char** argv,
                    AboutData const& ad,
                    po::options_description const& od )
    :
    Application( argc, argv, ad, od ),
    M_meshSize( 0.1 ),
    exporter( Exporter<mesh_type>::New( this->vm(), this->about().appName() ) )
{
    this->setOptions();
}
template<int Dim>
void
MyMesh<Dim>::setOptions()
{
    M_meshSize = this->vm()["hsize"].template as<double>();
}
template<int Dim>
typename MyMesh<Dim>::mesh_ptrtype
MyMesh<Dim>::createMesh()
{
    //# marker3 #
    mesh_ptrtype mesh( new mesh_type );
    //# endmarker3 #

    //# marker4 #
    GmshTensorizedDomain<convex_type::nDim,
                         convex_type::nOrder,
                         convex_type::nDim,
                         Simplex> geo;
    geo.setCharacteristicLength( M_meshSize );
    std::string fname = geo.generate( "mymesh" );
    //# endmarker4 #

    //# marker5 #
    ImporterGmsh<mesh_type> import( fname );
    mesh->accept( import );
    mesh->setComponents( MESH_PARTITION| MESH_UPDATE_FACES|MESH_UPDATE_EDGES);
    mesh->updateForUse();
    //# endmarker5 #

    return mesh;
}
template<int Dim>
void MyMesh<Dim>::run()
{
    if ( this->vm().count( "help" ) )
        {
            std::cout << this->optionsDescription() << "\n";
            return;
        }
    this->changeRepository( boost::format( "doc/tutorial/%1%/h_%2%/" )
                            % this->about().appName()
                            % M_meshSize
                            );

    mesh_ptrtype mesh = this->createMesh();

    //# marker6 #
    typedef bases<Lagrange<0,Scalar> > p0_basis_type;
    typedef FunctionSpace<mesh_type, p0_basis_type, Discontinuous> p0_space_type;
    boost::shared_ptr<p0_space_type> P0h( new p0_space_type( mesh ) );

    exporter->step(0)->setMesh( mesh );
    exporter->step(0)->add( "pid", regionProcess( P0h ) );

    exporter->save();
    //# endmarker6 #
}

//
// main function: entry point of the program
//
int main( int argc, char** argv )
{
    MyMesh<2> app( argc, argv, makeAbout(), makeOptions() );

    app.run();
}


