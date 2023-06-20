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
#define BOOST_TEST_MODULE test mmg gradedls

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
using namespace std::string_literals;

inline
po::options_description
makeOptions()
{
    // clang-format: off
    po::options_description desc_options( "test_Mmg options" );
    desc_options.add_options()
        ( "json", po::value<std::string>(), "json file" )
        ( "niter", po::value<int>()->default_value( 1 ), "niter" )
        ( "hsize", po::value<double>()->default_value( 0.1 ), "mesh size" )
        ( "remesh.verbose", po::value<int>()->default_value( 1 ), "mesh size" )
        ( "remesh.hmin", po::value<double>(), "Minimal mesh size " )
        ( "remesh.hmax", po::value<double>(), "Maximal mesh size " )
        ( "remesh.hsiz", po::value<double>(), "Constant mesh size " )
        ;
    // clang-format: on
    return desc_options.add( backend_options("I") ).add( backend_options("Iv") ).add( Feel::feel_options() );
}

/*_________________________________________________*
 * About
 *_________________________________________________*/

inline
AboutData
makeAbout()
{
    AboutData about( "test_remesh_graded" ,
                     "test_remesh_graded" ,
                     "0.1",
                     "test Mmg with graded levelset metric",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2020-2023 Universite de Strasbourg" );

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
    
    TestRemesh(std::string const& jname = "test_remesh_graded.json" )
    {
        
        std::string fname = Environment::expand(fmt::format("$cfgdir/{}",jname));
        BOOST_MESSAGE( fmt::format( "initialize test_remesh with json file: {}\n", fname) );
        if ( fs::exists(fs::path(fname) ) )
        {
            std::ifstream fin(fname.c_str());
            fin >> j_;
            BOOST_MESSAGE( fmt::format( "json: {}\n",j_.dump(1) ) );
        }
    }
    mesh_ptr_t setMesh( std::string const& fname = "", std::string const& e = "", std::vector<std::string> const& f = {}, std::vector<std::string> const& f_c = {} );
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
    nl::json j_;
    mesh_ptr_t mesh_;
    std::string required_elements_str_;
    std::vector<std::string> required_facets_str_;
    std::vector<std::string> check_facets_;
};
template <int Dim, int RDim>
typename TestRemesh<Dim, RDim>::mesh_ptr_t 
TestRemesh<Dim, RDim>::setMesh( std::string const& filename, std::string const& r_e, std::vector<std::string> const& r_f, std::vector<std::string> const& r_c )
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
    //std::vector<std::tuple<int,std::string>> iters{ {0,"0.2"}, {1,"0.1"}, {2,"0.05"} };
    std::vector<std::tuple<int, std::string>> iters{ { 0, "gradedls([RequiredBoundaryOfRequiredElements],0.1,0.25,0.25)" },
                                                     { 1, "gradedls([RequiredBoundaryOfRequiredElements],0.05,0.25,0.1)" },
                                                     { 2, "gradedls([RequiredBoundaryOfRequiredElements],0.025,0.25,0.1)" },
                                                     { 3, "gradedls([RequiredBoundaryOfRequiredElements],0.01,0.1,0.0.05)" } };
    for( auto [iter,metric]: iters) //ranges::views::ints( 0, niter ) )
    {
        if ( iter >= niter )
            break;
        auto Vh = createSpace<scalar,2>( mesh_ );
        auto Mh = createMetricSpace( mesh_ );
        auto met = Mh->element();
        bool created = createGradedLevelsetMetric( Mh, metric, 
                                                   [&met](auto const& met_ ){ met=met_; }, required_facets_str_ );
        if ( !created )
        {
            BOOST_MESSAGE( fmt::format("metric = {} not created",metric ) );
            throw std::invalid_argument( fmt::format("metric = {} not created",metric )  );
        }
        dump( mesh_, fmt::format( "dump mesh before remesh dim: {} iter: {}, metric: {}", Dim, iter, metric ) );
        auto [out,cpt] = remesh( _mesh=mesh_, _metric=metric, 
                                 _required_elts=std::vector<std::string>{{ required_elements_str_ }}, 
                                 _required_facets=required_facets_str_,
                                 _params=j_ );
        dump( out, fmt::format( "dump mesh before remesh dim: {} iter: {}, metric: {}", Dim, iter, metric ) );
        BOOST_CHECK_EQUAL(  nelements( markedelements(mesh_,required_elements_str_), true),nelements( markedelements(out,required_elements_str_),true) );
        BOOST_CHECK_EQUAL(  nelements( markedfaces(mesh_,required_facets_str_),true),nelements( markedfaces(out,required_facets_str_),true) );

        auto ein = exporter( _mesh = mesh_, _name = fmt::format( "initial-{}-{}", Dim, iter ) );
        auto eout = exporter( _mesh = out, _name = fmt::format("remeshed-{}-{}", Dim, iter) );
        eout->addRegions();
        ein->addRegions();
        ein->add( "metric", met );

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
BOOST_AUTO_TEST_CASE_TEMPLATE( testGradedMetric, T, dim_types )
{
    using namespace Feel;
    
    if ( Environment::isParallel() && T::first_type::value == 2 )
        return;
    double meshSize = doption("gmsh.hsize");

    if ( Environment::isMasterRank() )
    {
        BOOST_MESSAGE( "================================================================" );
        BOOST_MESSAGE( fmt::format( "[test_required_elements_and_facets] dim={}", T::first_type::value) );
    }
    TestRemesh<T::first_type::value, T::second_type::value> r("test_remesh_graded.json");
    r.setMesh( fmt::format("$cfgdir/domains_{}d.geo",T::first_type::value) );
    r.execute( ioption("niter") );
}

BOOST_AUTO_TEST_SUITE_END()
