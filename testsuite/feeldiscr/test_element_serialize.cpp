/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

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

#define BOOST_TEST_MODULE test_element_serialize
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelcore/serialization.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelvf/vf.hpp>

/** use Feel namespace */
using namespace Feel;
//using namespace Feel::vf;


inline
po::options_description
makeOptions()
{
    po::options_description testelementserializeoptions( "TestElementSerialize options" );
    testelementserializeoptions.add_options()
        ( "nb_element", po::value<int>()->default_value( 5 ), "number of elements of the vector" )
        ;
    return testelementserializeoptions.add( Feel::feel_options() );
}

template<int Dim,int Order>
class TestElementSerialize
{
public:

    //static const uint16_type Order = 1;

    typedef double value_type;
    typedef Simplex<Dim> convex_type;
    typedef Mesh<convex_type> mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;
    typedef bases<Lagrange<Order,Scalar> > basis_type;
    typedef FunctionSpace<mesh_type, basis_type> space_type;
    typedef typename space_type::element_type element_type;
    /**
     * Constructor
     */
    TestElementSerialize( std::string const& dbName, bool create_database )
        :
        M_dbName( dbName ),
        M_dbDirectory( fs::current_path()/dbName ),
        M_nb_element( 0 )
    {
        if ( create_database )
        {
            M_nb_element = ioption(_name="nb_element");
            double meshSize = 0.3;
            this->create_database( meshSize );
        }
        else
        {
            this->load_database();
        }
    }

    std::shared_ptr<space_type> functionSpace() const { return M_functionSpace; }
    element_type element() const { return M_element; }
    std::vector< std::shared_ptr<element_type> > const& vector_element() const { return M_vector_element; }

    void create_database( double meshSize );
    void load_database();

    fs::path dbPath() const { return M_dbDirectory/"database"; }



    friend class boost::serialization::access;

    BOOST_SERIALIZATION_SPLIT_MEMBER()

    template<class Archive>
    void save( Archive & ar, const unsigned int version ) const;
    template<class Archive>
    void load( Archive & ar, const unsigned int version ) ;

private:

    std::string M_dbName;
    fs::path M_dbDirectory;

    std::shared_ptr<space_type> M_functionSpace;
    element_type M_element;
    int M_nb_element;
    std::vector< std::shared_ptr<element_type> > M_vector_element;
};


template<int Dim,int Order>
void
TestElementSerialize<Dim,Order>::create_database( double meshSize )
{
    BOOST_TEST_MESSAGE( "DB will be created " );
    auto mesh = loadMesh(_mesh = new mesh_type,_h=meshSize );
    mesh->save( _name="mymesh",_path=M_dbDirectory,_type="text" );

    M_functionSpace = space_type::New( _mesh=mesh );
    M_element = M_functionSpace->element();
    M_vector_element.resize( M_nb_element );

    M_element.on( _range=elements(mesh), _expr=Px() );
    M_element.setName("element");
    for( int i = 0 ; i < M_nb_element ; i++ )
    {
        M_vector_element[i] = M_functionSpace->elementPtr();
        M_vector_element[i]->on( _range=elements(mesh), _expr=cst(i+3.14)*Px() );
        M_vector_element[i]->setName("element_"+std::to_string(i));
    }

    if ( !fs::exists( M_dbDirectory ) )
        fs::create_directories( M_dbDirectory );

    std::ofstream ofs ( this->dbPath() );
    CHECK( ofs ) << "can't write the db";
    boost::archive::text_oarchive oa( ofs );
    oa << *this;

}

template<int Dim,int Order>
void
TestElementSerialize<Dim,Order>::load_database()
{
    auto mesh = mesh_type::New();
    bool mesh_is_loaded = mesh->load( _name="mymesh",_path=M_dbDirectory,_type="text" );
    CHECK( mesh_is_loaded ) << "can't load the mesh";
    //this->load
    M_functionSpace = space_type::New( _mesh=mesh );
    M_element = M_functionSpace->element();


    std::ifstream ifs( this->dbPath() );
    CHECK( ifs ) << "can't read the db";
    boost::archive::text_iarchive ia( ifs );
    ia >> *this;
}


template< int Dim,int Order>
template<class Archive>
void
TestElementSerialize<Dim,Order>::save( Archive & ar, const unsigned int version ) const
{
    ar & BOOST_SERIALIZATION_NVP( M_nb_element );
    ar & BOOST_SERIALIZATION_NVP( M_element );

    for( int e=0; e<M_nb_element; e++ )
        ar & BOOST_SERIALIZATION_NVP( *(M_vector_element[e]) );

    //ar & BOOST_SERIALIZATION_NVP( M_vector_element );
}

template< int Dim,int Order>
template<class Archive>
void
TestElementSerialize<Dim,Order>::load( Archive & ar, const unsigned int version )
{
    ar & BOOST_SERIALIZATION_NVP( M_nb_element );

    ar & BOOST_SERIALIZATION_NVP( M_element );

    M_vector_element.resize( M_nb_element );

#if 1
    for(int e=0; e<M_nb_element; e++)
    {
        M_vector_element[e] = M_functionSpace->elementPtr();
        ar & BOOST_SERIALIZATION_NVP( *(M_vector_element[e]) );
    }
#else
    //problem when loading M_vector_element
    ar & BOOST_SERIALIZATION_NVP( M_vector_element );
#endif
}


template <typename DatabaseType>
void compareDatabases( DatabaseType const& db1, DatabaseType const& db2 )
{
    auto space1 = db1.functionSpace();
    auto space2 = db2.functionSpace();
    BOOST_CHECK( space1->nDof() > 0 );
    BOOST_CHECK( space1->nDof() == space2->nDof() );

    auto rangeElt = elements(space1->mesh());
    auto const& elt1 = db1.element();
    auto const& elt2 = db2.element();
    BOOST_CHECK( elt1.name() == elt2.name() );
    BOOST_CHECK_SMALL( normL2(_range=rangeElt,_expr=idv(elt1)-idv(elt2)), 1e-8 );

    auto const& vec_elt1 = db1.vector_element();
    auto const& vec_elt2 = db2.vector_element();
    BOOST_CHECK( vec_elt1.size() == vec_elt2.size() );
    for( int e = 0 ; e < vec_elt1.size(); e++ )
        BOOST_CHECK_SMALL( normL2(_range=rangeElt,_expr=idv(vec_elt1[e])-idv(vec_elt2[e])), 1e-8 );
}



FEELPP_ENVIRONMENT_WITH_OPTIONS( Feel::makeAboutDefault("test_element_serialize"), makeOptions() )

BOOST_AUTO_TEST_SUITE( element_serialize )

typedef boost::mpl::list<
    std::pair<boost::mpl::int_<1>,boost::mpl::int_<1>>,
    std::pair<boost::mpl::int_<1>,boost::mpl::int_<2>>,
    std::pair<boost::mpl::int_<2>,boost::mpl::int_<1>>,
    std::pair<boost::mpl::int_<2>,boost::mpl::int_<2>>,
    std::pair<boost::mpl::int_<3>,boost::mpl::int_<1>>,
    std::pair<boost::mpl::int_<3>,boost::mpl::int_<2>> > dim_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( MyElementSerializeCase, T, dim_types )
{
    static const int dim = T::first_type::value;
    static const int order = T::second_type::value;

    std::string dbName = (boost::format("db_%1%d_P%2%")%dim%order).str();
    TestElementSerialize<dim,order> dbc( dbName, true );
    TestElementSerialize<dim,order> dbl( dbName, false );
    compareDatabases( dbc, dbl );
}
BOOST_AUTO_TEST_SUITE_END()
