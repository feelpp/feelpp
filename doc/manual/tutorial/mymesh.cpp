/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

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
#include <feel/options.hpp>
#include <feel/feelcore/application.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feeltiming/tic.hpp>
#include <feel/feelfilters/gmsh.hpp>

#include <feel/feelfilters/exporter.hpp>

using namespace Feel;

inline
Feel::po::options_description
makeOptions()
{
    Feel::po::options_description mymeshoptions( "MyMesh options" );
    mymeshoptions.add_options()
    ( "hsize", Feel::po::value<double>()->default_value( 0.1 ), "mesh size" )
    ( "shape", Feel::po::value<std::string>()->default_value( "hypercube" ), "shape of the domain (either simplex or hypercube)" )
    ;

    // return the options mymeshoptions and the feel_options defined
    // internally by Feel
    return mymeshoptions.add( Feel::feel_options() );
}
inline
Feel::AboutData
makeAbout()
{
    Feel::AboutData about( "mymesh" ,
                           "mymesh" ,
                           "0.2",
                           "my first Feel application",
                           Feel::AboutData::License_GPL,
                           "Copyright (c) 2008-2010 Universite Joseph Fourier" );

    about.addAuthor( "Christophe Prud'homme",
                     "developer",
                     "christophe.prudhomme@ujf-grenoble.fr", "" );
    return about;
}

class A { public:
A():str(){}
A(std::string s):str(s){}
std::string str;
    friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
        {
            ar & str;
        }
};

template<int Dim>
class MyMesh: public Feel::Simget
{
public:

    //# marker1 #
    typedef Simplex<Dim> convex_type;
    //typedef Hypercube<Dim, 1,Dim> convex_type;
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

    Environment::changeRepository( boost::format( "doc/tutorial/%1%/%2%-%3%/h_%4%/" )
                                   % this->about().appName()
                                   % shape
                                   % Dim
                                   % meshSize );
    //Environment::setLogs( this->about().appName() );
    //# marker4 #
    auto mesh = mesh_type::New();

    tic();
    auto is_mesh_loaded = mesh->load( _name="mymesh",_path=".",_type="text" );
    toc( "load" );

    if ( !is_mesh_loaded )
    {
        tic();
        Log() << "Generating mesh using gmsh...\n";
        mesh = createGMSHMesh( _mesh=new mesh_type,
                               _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                               _desc=domain( _name=( boost::format( "%1%-%2%" ) % shape % Dim ).str() ,
                                             _shape=shape,
                                             _dim=Dim,
                                             _h=X[0] ),
                               _partitions=this->comm().size() );
        Log() << "Saving mesh...\n";
        mesh->save( _name="mymesh",_path=".",_type="text" );
        toc("generate+save");
    }
    else
    {
        tic();
        mesh->save( _name="mymesh2",_path=".",_type="text" );
        toc("save");
    }

    //# endmarker4 #

    int ne = std::distance( mesh->beginElementWithProcessId( this->comm().rank() ),
                            mesh->endElementWithProcessId( this->comm().rank() ) );
    Log() << "Local number of elements: " << ne << "\n";
    int gne;
    mpi::all_reduce( this->comm(), ne, gne, [] ( int x, int y )
    {
        return x + y;
    } );
    Log() << "Global number of elements: " << gne << "\n";

    //# marker62 #
    if ( exporter->doExport() )
    {
        Log() << "Exporting mesh\n";
        exporter->step( 0 )->setMesh( mesh );
        Log() << "Exporting regions\n";
        exporter->step( 0 )->addRegions();
        Log() << "Saving...\n";
        exporter->save();
    }
    //# endmarker62 #

#if 1

#if 1
    mpi::communicator world;
    Log() << "Sending/Receiving mesh...\n";
    if ( world.rank() == 0 )
    {
        boost::shared_ptr<A> str( new A("Hello") );
        world.send( 1, 10, mesh );
        //world.send( 1, 10, str);
        Log() << "done iSending\n";
    }
    else
    {
        Log() << "Receiving mesh from : " << 0 << "by " << world.rank() << "...\n";
        boost::shared_ptr<A> str;
        mesh_ptrtype mesh2;
        Log() << "Receiving mesh from : " << 0 << "by " << world.rank() << "...\n";
        Log() << "Receiving mesh from : " << 0 << "by " << world.rank() << "...\n";
        Log() << "Receiving mesh from : " << 0 << "by " << world.rank() << "...\n";
        //world.recv( 0, 10, str );
        world.recv( 0, 10, mesh2 );
        Log() << "Saving received mesh...\n";
        Log() << str->str << "...\n";
        //mesh2->save( _name="mymesh3", _type="text", _path="." );
    }
#else

    mpi::communicator world;

    if (world.rank() == 0) {
        world.send(1, 0, std::string("Hello"));
        std::string msg;
        world.recv(1, 1, msg);
        std::cout << msg << "!" << std::endl;
    } else {
        std::string msg;
        world.recv(0, 0, msg);
        std::cout << msg << ", ";
        std::cout.flush();
        world.send(0, 1, std::string("world"));
    }

#endif

#else
    if ( this->comm().size() > 1 )
    {
        Log() << "Sending meshes to proc 0\n";
        // send all meshes to proc 0
        if (  this->comm().rank() != 0 )
        {
            Log() << "gather mesh in 0\n";
            mpi::gather( this->comm(), mesh, 0 );
        }
        else
        {
            Log() << "gather mesh in 0\n";
            std::vector<mesh_ptrtype> mesh_vector;
            mpi::gather( this->comm(), mesh, mesh_vector, 0 );
            Log() << "mesh gathered in 0:" << mesh_vector.size() << "\n";
#if 0
            for( int i = 0; i < mesh_vector.size();++i )
            {
                exporter->step( i+1)->setMesh( mesh_vector[i] );
                exporter->save();
            }
#endif
        }

    }
#endif

}

//
// main function: entry point of the program
//
int main( int argc, char** argv )
{
    Feel::Environment env( argc, argv );

    Application app( argc, argv, makeAbout(), makeOptions() );

    if ( app.vm().count( "help" ) )
    {
        std::cout << app.optionsDescription() << "\n";
        return 0;
    }

    //app.add( new MyMesh<1>( app.vm(), app.about() ) );
    app.add( new MyMesh<2>( app.vm(), app.about() ) );
    //app.add( new MyMesh<3>( app.vm(), app.about() ) );

    app.run();
}


