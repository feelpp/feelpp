/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 16 Feb 2020

 Copyright (C) 2020 Feel++ Consortium

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
#define BOOST_TEST_MODULE test mmg

#include <range/v3/view/iota.hpp>
#include <feel/feelcore/testsuite.hpp>
#include <boost/hana/integral_constant.hpp>

#include <feel/feeldiscr/mesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feeldiscr/operatorinterpolation.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feelfilters/unitcube.hpp>
#include <feel/feeldiscr/createsubmesh.hpp>
#include <feel/feelfilters/exporter.hpp>
#include <feel/feelmesh/remesh.hpp>
#include <feel/feelpde/cg_laplacian.hpp>
#include <feel/feelmesh/dump.hpp>

using namespace Feel;
namespace hana =  boost::hana;
using namespace hana::literals;
inline
po::options_description
makeOptions()
{
    po::options_description desc_options( "test_Mmg options" );
    desc_options.add_options()
        ( "hsize", po::value<double>()->default_value( 0.1 ), "mesh size" )
        ( "remesh.verbose", po::value<int>()->default_value( 1 ), "mesh size" )
        ( "remesh.hmin", po::value<double>(), "Minimal mesh size " )
        ( "remesh.hmax", po::value<double>(), "Maximal mesh size " )
        ( "remesh.hsiz", po::value<double>(), "Constant mesh size " )
        
        ;
    return desc_options.add( backend_options("I") ).add( backend_options("Iv") ).add( Feel::feel_options() );
}

/*_________________________________________________*
 * About
 *_________________________________________________*/

inline
AboutData
makeAbout()
{
    AboutData about( "test_remesh" ,
                     "test_remesh" ,
                     "0.1",
                     "test Mmg",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2020 Universite de Strasbourg" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;
}


/**
 * @brief test remesh with various options
 * 
 * \tparam Dim real dimension
 * \tparam ScalarOrVectorial = 0 if scalar problem, 1 if vectorial
 * \tparam RDim real dimension
 */
template<int Dim, int RDim  = Dim>
class TestRemesh
{
    public:
    static constexpr int topo_dim = Dim;
    static constexpr int real_dim = RDim;
    static constexpr int scalar = 0;
    static constexpr int vectorial = 1;

    using mesh_t = Mesh<Simplex<Dim,1,RDim>>;
    using mesh_ptr_t = std::shared_ptr< mesh_t >;
    
    TestRemesh() = default;
    mesh_ptr_t setMesh( std::string const& fname = "", std::string const& e = "", std::string const& f = "", std::vector<std::string> const& f_c = {} );
    auto createMetricSpace( mesh_ptr_t m )
    {
        return Pch<1>( m );
    }
    template<int ScalarOrVectorial, int Order = 1>
    auto createSpace( mesh_ptr_t m ) 
    {
        if constexpr ( ScalarOrVectorial == scalar )
        {
            return Pch<Order>( m );
        }
        else if constexpr ( ScalarOrVectorial == vectorial )
        {
            return Pchv<Order>( m );
        }
    }
    std::vector<std::pair<std::string,std::string>>
    checkers()
    {
        if constexpr ( real_dim == 2 )
        {
            return std::vector<std::pair<std::string, std::string>>{ { "10", "{10,11}" }, { "0.3*x:x", "{0.3*x,0.4*x}:x" }, { "0.3*x+0.4*y:x:y", "{0.3*x+0.4*y,0.4*x+5*y}:x:y" }, { "0.3*x+0.4*y+x*x+2*y*y:x:y", "{0.3*x+0.4*y+x*x+2*y*y,0.4*x+5*y+x*x+2*y*y}:x:y" } };
        }
        else if constexpr ( real_dim == 3 )
        {
            //clang-format off
            return std::vector<std::pair<std::string, std::string>>{ { "10", "{10,11,12}" }, { "0.3*x:x", "{0.3*x,0.4*x,0.5*x}:x" }, { "0.3*x+0.4*y:x:y", "{0.3*x+0.4*y,0.4*x+5*y,0.4*x+5*y}:x:y" }, { "0.3*x+0.4*y:x:y", "{0.3*x+0.4*y,0.4*x+5*y,0.4*x+5*y}:x:y" }, { "0.3*x+0.4*y+10*z:x:y:z", "{0.3*x+0.4*y+0.5*z,0.4*x+5*y+6*z,0.4*x+5*y+7*z}:x:y:z" }, { "0.3*x+0.4*y+10*z+x*x+y*y-2*z*z:x:y:z", "{0.3*x+0.4*y+0.5*z+x*x+y*y-2*z*z,0.4*x+5*y+6*z+x*x+y*y-2*z*z,0.4*x+5*y+7*z+x*x+y*y-2*z*z}:x:y:z" } };
            //clang-format on
        }
    }
    int execute( int loop = 1 );
    mesh_ptr_t mesh_;
    std::string required_elements_str_;
    std::string required_facets_str_;
    std::vector<std::string> check_facets_;
};
template <int Dim, int RDim>
typename TestRemesh<Dim, RDim>::mesh_ptr_t 
TestRemesh<Dim, RDim>::setMesh( std::string const& filename, std::string const& r_e, std::string const& r_f, std::vector<std::string> const& r_c )
{
    auto create_mesh = [&]() {
        if constexpr ( topo_dim == 2 && real_dim == 2 )
        {
            if ( filename.empty() )
                return unitSquare();
            else
                return loadMesh( _mesh = new Mesh<Simplex<2>>, _filename = filename );
        }
        else if constexpr ( topo_dim == 2 && real_dim == 3 )
        {
            auto m = unitCube();
            auto mS = createSubmesh( _mesh = m, _range = boundaryfaces( m ) );
            return mS;
        }
        else if constexpr ( topo_dim == 3 )
        {
            if ( filename.empty() )
                return unitCube();
            else
                return loadMesh( _mesh = new Mesh<Simplex<3, 1>>( "mesh" ), _filename = filename );
        }
    };
    mesh_ = create_mesh();
    required_elements_str_=r_e;
    required_facets_str_ = r_f;
    check_facets_ = r_c;
    dump( mesh_ );
    return mesh_;
}
template <int Dim, int RDim>
int
TestRemesh<Dim, RDim>::execute( int niter )
{
    std::vector<std::tuple<int,std::string>> iters{ {0,"0.2"}, {1,"0.1"}, {2,"0.05"} };
    for( auto [iter,metric]: iters) //ranges::views::ints( 0, niter ) )
    {
        auto Vh = createSpace<scalar,2>( mesh_ );
        auto Wh = createSpace<vectorial,2>( mesh_ );
        
        auto Mh = createMetricSpace( mesh_ );
        auto met = Mh->element();
        //met.on( _range = elements( mesh_ ), _expr = expr( soption( "functions.s" ) ) );
        met.on( _range = elements( mesh_ ), _expr = expr( metric ) );

        auto r = remesher( mesh_, mesh_->markerName( required_elements_str_ ), mesh_->markerName( required_facets_str_ ) );

        r.setMetric( met );
        auto out = r.execute();
        out->updateForUse();
        dump( out );

        auto Mhr = createMetricSpace( out );
        auto Vhr = createSpace<scalar,2>( out );
        auto Whr = createSpace<vectorial,2>( out );
        auto interp = I( _domainSpace = Vh, _imageSpace = Vhr, _backend=backend(_name="I",_rebuild=true) );
        auto interpv = I( _domainSpace = Wh, _imageSpace = Whr, _backend=backend(_name="Iv",_rebuild=true) );

        auto ein = exporter( _mesh = mesh_, _name = fmt::format( "initial-{}", iter ) );
        auto eout = exporter( _mesh = out, _name = fmt::format("remeshed-{}", iter) );
        eout->addRegions();
        auto metout = Mhr->element( expr( metric ) );
        eout->add( "metricout", metout );
        ein->addRegions();
        ein->add( "metric", met );

        for( auto f : checkers() )
        {
            auto v = Vh->element( expr(f.first ) );
            ein->add( "v", v );
            auto iv = interp.operator()( v );
            auto vr = Vhr->element( expr( f.first ) );
            eout->add( "iv", iv );
            eout->add( "v", vr );
            double errv = normL2( _range = elements( out ), _expr = idv( vr ) - idv( iv ) );
            if ( Environment::isMasterRank() )
                BOOST_MESSAGE( fmt::format("L2 error norm {}: {}", f.first, errv ) );
            BOOST_CHECK_SMALL(errv, 1e-10);

            auto w = Wh->element( expr<real_dim, 1>( f.second ) );
            ein->add( "w", w );
            auto iw = interpv.operator()( w );
            auto wr = Whr->element( expr<real_dim,1>( f.second ) );
            eout->add( "iw", iw );
            eout->add( "w", wr );
            double errw = normL2( _range = elements( out ), _expr = idv( wr ) - idv( iw ) );
            if ( Environment::isMasterRank() )
                BOOST_MESSAGE( fmt::format( "L2 error norm {}: {}", f.second, errw ) );
            BOOST_CHECK_SMALL( errw, 1e-10 );
        }
        if ( !required_elements_str_.empty() )
        {
            double measinit = measure( _range = markedelements( Vh->mesh(), required_elements_str_ ) );
            double measremesh = measure( _range = markedelements( Vhr->mesh(), required_elements_str_ ) );
            if ( Environment::isMasterRank() )
                BOOST_MESSAGE( fmt::format( "check element set measure {} initial: {} remesh: {} error: {}", required_elements_str_, measinit, measremesh, measinit - measremesh ) );
            BOOST_CHECK_SMALL( measinit-measremesh, 1e-10 );
        }
        if ( !required_facets_str_.empty() )
        {
            double measinit = measure( _range = markedfaces( Vh->mesh(), required_facets_str_ ) );
            double measremesh = measure( _range = markedfaces( Vhr->mesh(), required_facets_str_ ) );
            if ( Environment::isMasterRank() )
                BOOST_MESSAGE( fmt::format( "check element set measure {} initial: {} remesh: {} error: {}", required_facets_str_, measinit, measremesh, measinit - measremesh ) );
            BOOST_CHECK_SMALL( measinit - measremesh, 1e-10 );
        }
        if ( !check_facets_.empty() )
        {
            for( auto m : check_facets_ )
            {
                double measinit = measure( _range = markedfaces( Vh->mesh(), m ) );
                double measremesh = measure( _range = markedfaces( Vhr->mesh(), m ) );
                if ( Environment::isMasterRank() )
                    BOOST_MESSAGE( fmt::format( "check element set measure {} initial: {} remesh: {} error: {}", m, measinit, measremesh, measinit - measremesh ) );
                BOOST_CHECK_SMALL( measinit - measremesh, 1e-10 );
            }
        }
        backend(_rebuild=true);
        //auto ur = cgLaplacianDirichlet( Vhr, std::tuple{ expr( "1" ), expr( "0" ), expr( "x+y:x:y" ) } );
        auto ur = cgLaplacianDirichlet( Vhr, std::tuple{ expr( "1" ), expr( "pi*pi*sin(pi*x):x" ), expr( "sin(pi*x):x" ) } );
        eout->add( "lap_u", ur );
        //eout->add( "lap_u_exact", expr( "x+y:x:y" ) );
        eout->add( "lap_u_exact", expr( "sin(pi*x):x" ) );
        //double errv = normL2( _range = elements( out ), _expr = idv( ur ) - expr( "x+y:x:y" ) );
        double errv = normL2( _range = elements( out ), _expr = idv( ur ) - expr( "sin(pi*x):x" ) );
        if ( Environment::isMasterRank() )
            BOOST_MESSAGE( fmt::format( "cglaplacian L2 error norm {}: {}", "x+y:x:y", errv ) );
        //BOOST_CHECK_SMALL( errv, 1e-10 );
        ein->save();
        eout->save();
        mesh_ = out;
    }
    return 1;
}
FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )


BOOST_AUTO_TEST_SUITE( mmg )

//typedef boost::mpl::list<boost::mpl::int_<2>,boost::mpl::int_<3> > dim_types;
typedef boost::mpl::list<
    std::pair<boost::mpl::int_<2>,boost::mpl::int_<2>>,
    //std::pair<boost::mpl::int_<2>,boost::mpl::int_<3>>,
    std::pair<boost::mpl::int_<3>,boost::mpl::int_<3>> > dim_types;
BOOST_AUTO_TEST_CASE_TEMPLATE( testnD, T, dim_types )
{
    using namespace Feel;
    
    if ( Environment::isParallel() && T::first_type::value == 2 )
        return;
    double meshSize = doption("gmsh.hsize");

    TestRemesh<T::first_type::value, T::second_type::value> r;
    r.setMesh();
    r.execute(3);

}

typedef boost::mpl::list<
    std::pair<boost::mpl::int_<2>, boost::mpl::int_<2>>,
    std::pair<boost::mpl::int_<3>, boost::mpl::int_<3>>>
    dim_types;
BOOST_AUTO_TEST_CASE_TEMPLATE( testnDMat, T, dim_types )
{
    using namespace Feel;

    if ( Environment::isParallel() && T::first_type::value == 2 )
        return;
    TestRemesh<T::first_type::value, T::second_type::value> r;
    if constexpr ( T::first_type::value == 2 )
    {
        r.setMesh( "$datadir/geo/square_squareHole.geo", "","BdryFixedDiscr" );
    }
    else
    {
        r.setMesh( "$datadir/geo/cube_cubeHole.geo", "","BdryFixedDiscr" );
    }
    r.execute( 1 );
}
typedef boost::mpl::list<
    std::pair<boost::mpl::int_<2>, boost::mpl::int_<2>>,
    std::pair<boost::mpl::int_<3>, boost::mpl::int_<3>>>
    dim_types;
BOOST_AUTO_TEST_CASE_TEMPLATE( testnDMatTwo, T, dim_types )
{
    using namespace Feel;

    if ( Environment::isParallel() && T::first_type::value == 2 )
        return;
    TestRemesh<T::first_type::value, T::second_type::value> r;
    if constexpr ( T::first_type::value == 2 )
    {
        r.setMesh( "$datadir/data/geo/square_twoMaterials.geo", "MatTwo", "BdryFixedDiscr" );
    }
    else
    {
        // r.setMesh( "$datadir/geo/cube_twoMaterials.geo", "MatTwo", "BdryFixedDiscr" );
        //r.setMesh( "$datadir/geo/cube_twoMaterials.geo", "MatTwo", "", { "BdryFixedDiscr"} );
        r.setMesh( "$datadir/geo/cube_twoMaterials.geo", "MatTwo");
    }
    r.execute( 1 );
}

BOOST_AUTO_TEST_SUITE_END()
