/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   Date: 2011-07-09

   Copyright (C) 2011 Universite Joseph Fourier (Grenoble I)

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
#define USE_BOOST_TEST 1

#define BOOST_TEST_MODULE integration_opt testsuite

#if defined(USE_BOOST_TEST)
#include <boost/test/unit_test.hpp>
using boost::unit_test::test_suite;
#endif

#include <feel/options.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feelfilters/gmsh.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelfilters/geotool.hpp>


using namespace Feel;
using namespace Feel::vf;

#if USE_BOOST_TEST

BOOST_AUTO_TEST_SUITE( integration_opt )
Feel::Environment env( boost::unit_test::framework::master_test_suite().argc,
                       boost::unit_test::framework::master_test_suite().argv );

typedef boost::mpl::list<boost::mpl::pair<mpl::int_<2>,mpl::int_<2> >,
                         boost::mpl::pair<mpl::int_<2>,mpl::int_<4> >,
                         boost::mpl::pair<mpl::int_<3>,mpl::int_<2> >,
                         boost::mpl::pair<mpl::int_<3>,mpl::int_<4> >
                         > dim_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( integration_opt, T, dim_types )
{
    double hsize = 2;
    // Declare the supported options.
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("hsize", po::value<double>(&hsize)->default_value( 2 ), "h size")
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(boost::unit_test::framework::master_test_suite().argc,
                                     boost::unit_test::framework::master_test_suite().argv,
                                     desc), vm);
    po::notify(vm);

    using namespace Feel;
    using namespace Feel::vf;
    typedef Mesh<Simplex<T::first::value,T::second::value> > mesh_type;

    auto mesh = createGMSHMesh( _mesh=new mesh_type,
                                _desc=domain( _name=(boost::format( "ellipsoid-%1%-%2%" ) % T::first::value % T::second::value).str() ,
                                              _usenames=true,
                                              _shape="ellipsoid",
                                              _dim=T::first::value,
                                              _order=T::second::value,
                                              _h=hsize ),
                                _update=MESH_CHECK|MESH_UPDATE_FACES|MESH_UPDATE_EDGES );
    //straightenMesh( _mesh=mesh );

    //std::cout << "read mesh\n" << std::endl;
    auto i1 = integrate( _range=elements(mesh), _quad=_Q<5>(), _expr=cst(1.), _geomap=GEOMAP_HO ).evaluate()( 0, 0 );
    auto i2 = integrate( _range=elements(mesh), _quad=_Q<5>(), _expr=cst(1.), _geomap=GEOMAP_OPT ).evaluate()( 0, 0 );
    BOOST_CHECK_CLOSE( i1, i2, 1e-12 );
    BOOST_TEST_MESSAGE( "ho = " << std::scientific << std::setprecision( 16 ) << i1 << "\n" <<
                        "opt = " << std::scientific << std::setprecision( 16 ) << i2 << "\n" << "\n");

}
BOOST_AUTO_TEST_SUITE_END()

#else

int
main( int argc, char** argv )
{

}

#endif
