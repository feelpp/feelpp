/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Stephane Veys <stephane.veys@imag.fr>
       Date: 2012-06-30

  Copyright (C) 2008-2010 Universite Joseph Fourier (Grenoble I)

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
   \file test_element_serialize.cpp
   \author Stephane Veys <stephane.veys@imag.fr>
   \date 2012-06-30
*/

#define BOOST_TEST_MODULE test_element_serialize
#include <testsuite/testsuite.hpp>

#include <fstream>

#include <boost/tuple/tuple.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/split_member.hpp>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <feel/feelcore/serialization.hpp>
#include <feel/feeldiscr/region.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/options.hpp>
#include <feel/feeltiming/tic.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelpoly/polynomialset.hpp>
#include <feel/feelvf/vf.hpp>

/** use Feel namespace */
using namespace Feel;
using namespace Feel::vf;


inline
po::options_description
makeOptions()
{
    po::options_description testelementserializeoptions( "TestElementSerialize options" );
    testelementserializeoptions.add_options()
        ( "hsize", po::value<double>()->default_value( 0.2 ), "mesh size" )
        ( "shape", Feel::po::value<std::string>()->default_value( "simplex" ), "shape of the domain (either simplex or hypercube)" )
        ( "nb_element", po::value<int>()->default_value( 2 ), "number of elements of the vector" )
        ;
    return testelementserializeoptions.add( Feel::feel_options() );
}


inline
AboutData
makeAbout()
{
    AboutData about( "test_element_serialize" ,
                     "test_element_serialize" ,
                     "0.2",
                     "nD(n=1,2,3) test serialization of element_type ( single and vector )",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2008-2009 Universite Joseph Fourier" );

    about.addAuthor( "Stephane Veys", "developer", "stephane.veys@imag.fr", "" );
    return about;

}


template<int Dim>
class TestElementSerialize
    :
public Simget
{
    typedef Simget super;
public:

    static const uint16_type Order = 1;

    typedef double value_type;
    typedef Simplex<Dim> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;
    typedef bases<Lagrange<Order,Scalar> > basis_type;
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef typename space_type::element_type element_type;
    /**
     * Constructor
     */
    TestElementSerialize( bool rebuild_database )
        :
        super(),
        M_meshSize( this->vm()["hsize"].template as<double>() ),
        M_shape( this->vm()["shape"].template as<std::string>() ),
        M_nb_element( this->vm()["nb_element"].template as<int>() ),
        M_rebuild_database( rebuild_database )
    {
    }

    void run();
    void run( const double* X, unsigned long P, double* Y, unsigned long N );

    void saveDB();
    void loadDB();
    fs::path dbPath();

    friend class boost::serialization::access;

    BOOST_SERIALIZATION_SPLIT_MEMBER()

    template<class Archive>
    void save( Archive & ar, const unsigned int version ) const;
    template<class Archive>
    void load( Archive & ar, const unsigned int version ) ;
    bool existDB();
    void setRebuildDatabase( bool b);

private:

    double M_meshSize;
    std::string M_shape;
    element_type M_element;
    element_type M_element_temp;
    int M_nb_element;
    bool M_rebuild_database;
    std::vector< element_type > M_vector_element;
    std::vector< element_type > M_rebuilt_vector_element;
}; // TestElementSerialize


template<int Dim>
void
TestElementSerialize<Dim>::run()
{
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Execute TestElementSerialize<" << Dim << ">\n";
    std::vector<double> X( 2 );
    X[0] = M_meshSize;

    if ( M_shape == "hypercube" )
        X[1] = 1;

    else // default is simplex
        X[1] = 0;

    std::vector<double> Y( 3 );
    run( X.data(), X.size(), Y.data(), Y.size() );
}

template<int Dim>
void
TestElementSerialize<Dim>::run( const double* X, unsigned long P, double* Y, unsigned long N )
{
    if ( X[1] == 0 ) M_shape = "simplex";

    if ( X[1] == 1 ) M_shape = "hypercube";

#if 1
    if ( !this->vm().count( "nochdir" ) )

        Environment::changeRepository( boost::format( "testsuite/feeldiscr/%1%/%2%-%3%/h_%4%/" )
                                       % this->about().appName()
                                       % M_shape
                                       % Dim
                                       % M_meshSize );
#endif


    auto mesh = mesh_type::New();

    auto is_mesh_loaded = mesh->load( _name="mymesh",_path=".",_type="text" );

    if ( ! is_mesh_loaded )
    {
        mesh = createGMSHMesh( _mesh=new mesh_type,
                               _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                               _desc=domain( _name=( boost::format( "%1%-%2%-%3%" ) % M_shape % Dim % 1 ).str() ,
                                             _usenames=true,
                                             _shape=M_shape,
                                             _dim=Dim,
                                             _h=X[0],
                                             _xmin=0.,
                                             _xmax=1.,
                                             _ymin=0.,
                                             _ymax=1.,
                                             _zmin=0.,
                                             _zmax=1. ) );
        mesh->save( _name="mymesh",_path=".",_type="text" );
    }


    if( ! existDB() || M_rebuild_database )
    {
        auto Xh = space_type::New( mesh );
        M_element = Xh->element();
        M_vector_element.resize( M_nb_element );
        std::cout<<"DB will be created "<<std::endl;
        M_element = vf::project( Xh , elements(mesh), cst( 1 ) );
        M_element.setName("element");
        for( int i = 0 ; i < M_nb_element ; i++ )
        {
            M_vector_element[i] = Xh->element();
            M_vector_element[i] = M_element ;
            M_vector_element[i].setName("element_"+(boost::format("%1%") %i ).str());
           std::cout<<"save element "<<M_vector_element[i].name()<<std::endl;
        }
        M_element_temp = Xh->element();
        this->saveDB();
        std::cout<<"DB successfully created "<<std::endl;
    }
    else
    {
        std::cout<<"DB will be loaded"<<std::endl;
        this->loadDB();
        std::cout<<"DB successfully loaded"<<std::endl;


        double max = M_element.max();
        double min = M_element.min();
        BOOST_CHECK_CLOSE( max, 1.0, 2e-1 );
        BOOST_CHECK_CLOSE( min, 1.0, 2e-1 );

        for( int e = 0 ; e < M_nb_element; e++ )
        {
            max = M_rebuilt_vector_element[e].max();
            min = M_rebuilt_vector_element[e].min();
            BOOST_CHECK_CLOSE( max, 1.0, 2e-1 );
            BOOST_CHECK_CLOSE( min, 1.0, 2e-1 );
        }

    }




} // TestElementSerialize::run

template< int Dim>
bool
TestElementSerialize<Dim>::existDB()
{
    return fs::exists( this->dbPath() / "database" );
}

template< int Dim>
template<class Archive>
void
TestElementSerialize<Dim>::save( Archive & ar, const unsigned int version ) const
{
    ar & BOOST_SERIALIZATION_NVP( M_nb_element );
    ar & BOOST_SERIALIZATION_NVP( M_element );

    for( int e=0; e<M_nb_element; e++ )
        ar & BOOST_SERIALIZATION_NVP( M_vector_element[e] );

    //ar & BOOST_SERIALIZATION_NVP( M_vector_element );
}

template< int Dim>
template<class Archive>
void
TestElementSerialize<Dim>::load( Archive & ar, const unsigned int version )
{

    auto mesh = mesh_type::New();
    auto is_mesh_loaded = mesh->load( _name="mymesh",_path=".",_type="text" );

    auto Xh = space_type::New( mesh );
    M_element = Xh->element();

    ar & BOOST_SERIALIZATION_NVP( M_nb_element );

    ar & BOOST_SERIALIZATION_NVP( M_element );

    M_rebuilt_vector_element.resize( M_nb_element );
    M_vector_element.resize( M_nb_element );

#if 1
    for(int e=0; e<M_nb_element; e++)
    {
        M_element_temp = Xh->element();
        M_element_temp.setName("element_"+(boost::format("%1%") %e ).str() );
        ar & BOOST_SERIALIZATION_NVP( M_element_temp );
        M_rebuilt_vector_element[e] = M_element_temp;

        M_vector_element[e] = Xh->element();
    }
#else
    //problem when loading M_vector_element
    ar & BOOST_SERIALIZATION_NVP( M_vector_element );
#endif
}


template< int Dim>
fs::path
TestElementSerialize<Dim>::dbPath()
{
    std::string localpath = ( boost::format( "%1%/testsuite/feeldiscr/%2%/%3%-%4%/h_%5%/db/" )
                              % Feel::Environment::rootRepository()
                              % this->about().appName()
                              % M_shape
                              % Dim
                              % M_meshSize).str();

    fs::path rep_path = localpath;

    return rep_path ;
}


template< int Dim >
void
TestElementSerialize<Dim>::saveDB()
{

    fs::path path = this->dbPath();
    fs::create_directories( path );

    fs::ofstream ofs ( this->dbPath() / "database" );

    if ( ofs )
    {
        boost::archive::text_oarchive oa( ofs );
        oa << *this;
    }
    else
        throw std::logic_error(" can't use ofs in saveDB() " );
}

template< int Dim >
void
TestElementSerialize<Dim>::loadDB()
{

    fs::path db = this->dbPath();

    fs::ifstream ifs( db / "database" );

    if ( ifs )
    {
        boost::archive::text_iarchive ia( ifs );
        ia >> *this;
    }
}


template < int Dim >
void
TestElementSerialize<Dim>::setRebuildDatabase( bool b )
{
    M_rebuild_database = b;
}

/**
 * main code
 */

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );

BOOST_AUTO_TEST_SUITE( element_serialize )

BOOST_AUTO_TEST_CASE( MyElementSerializeCase )
{
    Application app;

    //the first one :  create database ( if doesn't exist )
    app.add( new TestElementSerialize<1>( true ) );
    //the second one : load database
    app.add( new TestElementSerialize<1>(false ) );

    app.add( new TestElementSerialize<2>(true ) );
    app.add( new TestElementSerialize<2>(false ) );

    app.add( new TestElementSerialize<3>(true ) );
    app.add( new TestElementSerialize<3>(false ) );

    app.run();

}
BOOST_AUTO_TEST_SUITE_END()


