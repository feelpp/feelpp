/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 17 mai 2016
 
 Copyright (C) 2016 Feel++ Consortium
 
 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.
 
 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */
#include <sstream>

#include <boost/timer.hpp>
// Boost.Test
// make sure that the init_unit_test function is defined by UTF
//#define BOOST_TEST_MAIN
// give a name to the testsuite
#define BOOST_TEST_MODULE mesh filter testsuite
#include <testsuite.hpp>
#include <boost/mpl/list.hpp>

#include <feel/feelcore/feel.hpp>


using boost::unit_test::test_suite;

#include <feel/feelfilters/exporter.hpp>
#include <feel/feeldiscr/mesh.hpp>
#include <feel/feelfilters/loadmesh.hpp>



using namespace Feel;

inline
po::options_description
makeOptions()
{
    po::options_description desc_options( "test_importerarm ACUSIM Raw Mesh importer options" );
    desc_options.add_options()
        ( "input.mesh.filename", po::value<std::string>()->default_value(""), "input.mesh.filename" )

        ;
    return desc_options.add( Feel::feel_options() );
}

/*_________________________________________________*
 * About
 *_________________________________________________*/

inline
AboutData
makeAbout()
{
    AboutData about( "test_importerarm" ,
                     "test_importerarm" ,
                     "0.1",
                     "test_importerarm",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2016 Feel++ Consortium" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;
}


template<int Dim, template <uint16_type,uint16_type,uint16_type> class Entity = Simplex>
void
checkCreateGmshMesh( std::string const& shape, std::string const& convex = "Simplex" )
{
    typedef Mesh<Entity<Dim,1,Dim> > mesh_type;
    typedef std::shared_ptr<mesh_type> mesh_ptrtype;

    typedef Mesh<Simplex<3>> mesh_type;
    //size_type updateComponentsMesh = MESH_UPDATE_ELEMENTS_ADJACENCY|MESH_NO_UPDATE_MEASURES;
    size_type updateComponentsMesh = MESH_UPDATE_FACES_MINIMAL;
    //size_type updateComponentsMesh = MESH_UPDATE_FACES_MINIMAL|MESH_UPDATE_EDGES;
    tic();


    std::string msh = soption(_name="arm.filename");

    auto mesh = loadMesh(_mesh=new mesh_type,
                         _filename=msh,
                         _update=updateComponentsMesh );
    if ( Environment::isMasterRank() )
        std::cout << "loadMesh done" << std::endl;
    toc("loadMesh",FLAGS_v>0);

    BOOST_CHECK_NE( nelements(elements(mesh),true), 0 );
}

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )


BOOST_AUTO_TEST_SUITE( armsuite )


    //typedef boost::mpl::list<boost::mpl::int_<1>,boost::mpl::int_<2>,boost::mpl::int_<3> > dim_types;
typedef boost::mpl::list<boost::mpl::int_<3> > dim_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( gmshsimplex, T, dim_types )
{
    checkCreateGmshMesh<T::value>( "simplex" );
}



BOOST_AUTO_TEST_SUITE_END()
