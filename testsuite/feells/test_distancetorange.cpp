/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Thibaut Metivet <thibaut.metivet@inria.fr>
 Date: 13 Sept 2021

 Copyright (C) 2021 Feel++ Consortium

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
#define BOOST_TEST_MODULE test distance

#include <feel/feelcore/testsuite.hpp>
#include <boost/hana/integral_constant.hpp>
#include <fmt/core.h>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelfilters/unitcube.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feells/distancetorange.hpp>

using namespace Feel;
namespace hana =  boost::hana;
using namespace hana::literals;

inline
po::options_description
makeOptions()
{
    po::options_description desc_options( "test_distancetorange options" );
    desc_options.add_options()
        ( "hsize", po::value<double>()->default_value( 0.1 ), "mesh size" )
        ( "remesh.verbose", po::value<int>()->default_value( 1 ), "mesh size" )
        ( "remesh.hmin", po::value<double>(), "Minimal mesh size " )
        ( "remesh.hmax", po::value<double>(), "Maximal mesh size " )
        ( "remesh.hsiz", po::value<double>(), "Constant mesh size " )
        
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
    AboutData about( "test_distancetorange" ,
                     "test_distancetorange" ,
                     "0.1",
                     "test DistanceToRange",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2021 Feel++ Consortium" );

    about.addAuthor( "Thibaut Metivet", "developer", "thibaut.metivet@inria.fr", "" );
    return about;
}


/**
 * @brief test DistanceToRange with various options
 * 
 * \tparam Dim dimension
 */
template<uint16_type Dim, uint16_type Order>
class TestDistanceToRange
{
    public:
        using mesh_type = Mesh<Simplex<Dim>>;
        using mesh_ptrtype = std::shared_ptr< mesh_type >;
        
        TestDistanceToRange() = default;

        mesh_ptrtype setMesh( std::string const& fname = "" );
        mesh_ptrtype const& mesh() const { return M_mesh; }

        int run();

    private:
        mesh_ptrtype M_mesh;
};

template<uint16_type Dim, uint16_type Order>
typename TestDistanceToRange<Dim, Order>::mesh_ptrtype 
TestDistanceToRange<Dim, Order>::setMesh( std::string const& filename )
{
    auto create_mesh = [&]() {
        if constexpr ( Dim == 2 )
        {
            if ( filename.empty() )
                return unitSquare();
            else
                return loadMesh( _mesh = new Mesh<Simplex<2>>( "mesh" ), _filename = filename );
        }
        else if constexpr ( Dim == 3 )
        {
            if ( filename.empty() )
                return unitCube();
            else
                return loadMesh( _mesh = new Mesh<Simplex<3, 1>>( "mesh" ), _filename = filename );
        }
    };
    M_mesh = create_mesh();
    return M_mesh;
}

template<uint16_type Dim, uint16_type Order>
int
TestDistanceToRange<Dim, Order>::run()
{
    if( !M_mesh )
        this->setMesh();
    auto const& mesh = this->mesh();
    auto Vh = Pch<Order>( mesh );

    // Exporter
    auto exp = exporter( _mesh = mesh, _name = fmt::format( "distance_{}d_o{}", Dim, Order ) );
    exp->addRegions();

    // Distance to boundary
    auto distToBoundary = distanceToRange( _space=Vh, _range=boundaryfaces( mesh ) );
    exp->add( "distToBoundary", distToBoundary );

    // Narrow-band distance to boundary 
    auto distToBoundaryNarrowBand = distanceToRange( 
            _space=Vh, _range=boundaryfaces( mesh ), 
            _max_distance=3.*mesh->hAverage() );
    exp->add( "distToBoundaryNarrowBand", distToBoundaryNarrowBand );

    // Distance to boundary with narrow-band and stride
    auto distToBoundaryNarrowBandStride = distanceToRange( 
            _space=Vh, _range=boundaryfaces( mesh ), 
            _max_distance=3.*mesh->hAverage(),
            _stride=2.*mesh->hAverage()
            );
    exp->add( "distToBoundaryNarrowBandStride", distToBoundaryNarrowBandStride );

    // Export
    exp->save();

    return 1;
}

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )


BOOST_AUTO_TEST_SUITE( distance )

typedef boost::mpl::list< boost::mpl::int_<2>, boost::mpl::int_<3> > dim_types;
//typedef boost::mpl::list< boost::mpl::int_<2> > dim_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( testnD, T, dim_types )
{
    using namespace Feel;

    TestDistanceToRange<T::value, 1> test;
    test.setMesh();
    test.run();

}

BOOST_AUTO_TEST_SUITE_END()
