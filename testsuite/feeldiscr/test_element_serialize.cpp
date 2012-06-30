/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Abdoulaye Samake <abdoulaye.samake@e.ujf-grenoble.fr>
       Date: 2011-08-10

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
   \file test.cpp
   \author Stephane Veys <stephane.veys@imag.fr>
   \date 2011-08-10
*/

#define BOOST_TEST_MODULE test_element_serialize
#include <boost/test/unit_test.hpp>
using boost::unit_test::test_suite;

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
#include <feel/feelfilters/gmsh.hpp>
#include <feel/options.hpp>
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
                     "nD(n=1,2,3) TestLift on simplices or simplex products",
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
    typedef bases<Lagrange<Order,Scalar> > basis_type;
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef typename space_type::element_type element_type;

    /**
     * Constructor
     */
    TestElementSerialize( po::variables_map const& vm, AboutData const& about )
        :
        super( vm, about ),
        meshSize( this->vm()["hsize"].template as<double>() ),
        shape( this->vm()["shape"].template as<std::string>() )
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

private:

    double meshSize;
    std::string shape;
    element_type element;
}; // TestElementSerialize


template<int Dim>
void
TestElementSerialize<Dim>::run()
{
    std::cout << "------------------------------------------------------------\n";
    std::cout << "Execute TestElementSerialize<" << Dim << ">\n";
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
TestElementSerialize<Dim>::run( const double* X, unsigned long P, double* Y, unsigned long N )
{
    if ( X[1] == 0 ) shape = "simplex";

    if ( X[1] == 1 ) shape = "hypercube";

    //if ( !this->vm().count( "nochdir" ) )

        Environment::changeRepository( boost::format( "testsuite/feeldiscr/%1%/%2%-%3%/h_%4%/" )
                                       % this->about().appName()
                                       % shape
                                       % Dim
                                       % meshSize );

    auto mesh = createGMSHMesh( _mesh=new mesh_type,
                                _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES|MESH_RENUMBER,
                                _desc=domain( _name=( boost::format( "%1%-%2%-%3%" ) % shape % Dim % 1 ).str() ,
                                        _usenames=true,
                                        _shape=shape,
                                        _dim=Dim,
                                        _h=X[0],
                                        _xmin=0.,
                                        _xmax=1.,
                                        _ymin=0.,
                                        _ymax=1.,
                                        _zmin=0.,
                                        _zmax=1. ) );



    auto Xh = space_type::New( mesh );
    element = Xh->element();

    if( ! existDB() )
    {
        std::cout<<"DB will be created "<<std::endl;
        element.setName("element");
        element = vf::project( Xh , elements(mesh), cst( 1 ) );
        this->saveDB();
        std::cout<<"DB successfully created "<<std::endl;
    }
    else
    {
        std::cout<<"DB will be loaded"<<std::endl;
        this->loadDB();
        std::cout<<"DB successfully loaded"<<std::endl;
    }

    double max = element.max();
    double min = element.min();
    BOOST_CHECK_CLOSE( max, 1.0, 2e-1 );
    BOOST_CHECK_CLOSE( min, 1.0, 2e-1 );


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
    ar & BOOST_SERIALIZATION_NVP( element );
}

template< int Dim>
template<class Archive>
void
TestElementSerialize<Dim>::load( Archive & ar, const unsigned int version )
{
    ar & BOOST_SERIALIZATION_NVP( element );
}


template< int Dim>
fs::path
TestElementSerialize<Dim>::dbPath()
{
    std::string localpath = ( boost::format( "%1%/testsuite/feeldiscr/%2%/%3%-%4%/h_%5%/db/" )
                              % Feel::Environment::rootRepository()
                              % this->about().appName()
                              % shape
                              % Dim
                              % meshSize).str();

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


/**
 * main code
 */


BOOST_AUTO_TEST_SUITE( element_serialize )
Environment env( boost::unit_test::framework::master_test_suite().argc,
                 boost::unit_test::framework::master_test_suite().argv );
BOOST_AUTO_TEST_CASE( MyElementSerializeCase )
{

    Feel::Environment env( boost::unit_test::framework::master_test_suite().argc,
                           boost::unit_test::framework::master_test_suite().argv );

    Application app( boost::unit_test::framework::master_test_suite().argc,
                     boost::unit_test::framework::master_test_suite().argv, makeAbout(), makeOptions() );

    if ( app.vm().count( "help" ) )
        std::cout << app.optionsDescription() << "\n";

    //the first one :  create database ( if doesn't exist )
    app.add( new TestElementSerialize<1>( app.vm(), app.about() ) );
    //the second one : load database
    app.add( new TestElementSerialize<1>( app.vm(), app.about() ) );

    app.add( new TestElementSerialize<2>( app.vm(), app.about() ) );
    app.add( new TestElementSerialize<2>( app.vm(), app.about() ) );

    app.add( new TestElementSerialize<3>( app.vm(), app.about() ) );
    app.add( new TestElementSerialize<3>( app.vm(), app.about() ) );

    app.run();

}
BOOST_AUTO_TEST_SUITE_END()


