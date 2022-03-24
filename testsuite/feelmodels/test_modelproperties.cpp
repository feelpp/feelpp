/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 15 Mar 2015

 Copyright (C) 2015 Feel++ Consortium

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
#define BOOST_TEST_MODULE model properties testsuite

#include <feel/feelcore/testsuite.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feelmodels/modelproperties.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/vf.hpp>

using namespace Feel;

inline
po::options_description
makeOptions()
{
    po::options_description modelopt( "test model properties options" );
    modelopt.add_options()
        ("json_filename" , Feel::po::value<std::string>()->default_value( "$cfgdir/test.feelpp" ),
         "json file" )
        ;
    return  modelopt.add( Feel::feel_options() ) ;
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_modelproperties",
                     "test_modelproperties" );
    return about;

}

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )

BOOST_AUTO_TEST_SUITE( modelproperties )

BOOST_AUTO_TEST_CASE( test_materials )
{
    auto mesh = loadMesh(_mesh=new Mesh<Simplex<3> >);
    auto Xh = Pch<1>(mesh);
    auto g = Xh->element();
    auto d = Xh->element();

    ModelProperties model_props;
    model_props.setupFromFilenameOption();

    auto mats = model_props.materials();
    for ( auto const& matPair : mats )
    {
        auto mat = matPair.second;
        auto physics = mat.physics();
        auto name = mat.name();
        auto meshMarkers = mat.meshMarkers();

        auto rho = mat.property("rho");
        auto eta = mat.property("eta");
        auto f = mat.property("f");
#if 0
        auto fParam = mat.getScalar("f",{{"g",1.}});
        auto hList = mat.getScalar("h",{"g","t"},{idv(g),idv(d)});
        auto fVec = mat.getScalar("f",std::vector<std::string>(1,"g"), std::vector<decltype(idv(g))>(1,idv(g)));
        auto hParams = mat.getScalar("h",{"g"},{idv(g)}, {{"t",1.}});

        auto nu = mat.getVector<3>( "nu" );
        auto curlnu = curl(nu);
        auto muPair = mat.getVector<2>("mu",{"t",1.});
        auto muMap = mat.getVector<2>("mu", {{"t",1.}});

        auto chi = mat.getMatrix<2>( "chi" );
        auto xhiPair = mat.getMatrix<3>( "xhi", {"t",2.} );
        auto xhiMap = mat.getMatrix<3,3>( "xhi", {{"t",3.}});
#endif
#if 0
        Feel::cout << "properties for " << matPair.first << std::endl;
        Feel::cout << "\t" << name << std::endl;
        Feel::cout << "\thas " << physics.size() << " physics:" << std::endl;
        for( auto const& p : physics )
            Feel::cout << "\t\t" << p << std::endl;
        Feel::cout << "\thas " << meshMarkers.size() << " markers:" << std::endl;
        for( auto const& p : meshMarkers )
            Feel::cout << "\t\t" << p << std::endl;

        Feel::cout << "\t" << rhoInt << std::endl;
        Feel::cout << "\t" << etaDouble << std::endl;
        Feel::cout << "\t" << rho << std::endl;
        Feel::cout << "\t" << nu << std::endl;
        Feel::cout << "\t" << curlnu << std::endl;
        Feel::cout << "\t" << muPair << std::endl;
        Feel::cout << "\t" << muMap << std::endl;
        Feel::cout << "\t" << chi << std::endl;
        Feel::cout << "\t" << xhiPair << std::endl;
        Feel::cout << "\t" << xhiMap << std::endl;
#endif
    }
#if 0
    Feel::cout << mats.materialWithPhysic("electro").size() << " materials with electro physic" << std::endl;
    Feel::cout << mats.materialWithPhysic("thermo").size() << " materials with thermo physic" << std::endl;
#endif
    BOOST_CHECK_EQUAL(mats.materialWithPhysic("electro").size(), 2);
    BOOST_CHECK_EQUAL(mats.materialWithPhysic("thermo").size(), 1);

}

BOOST_AUTO_TEST_CASE( test_parameters )
{
    ModelProperties model_props;
    model_props.setupFromFilenameOption();

    auto param = model_props.parameters();
    BOOST_CHECK_EQUAL( param.size(), 3);
    for ( auto const& pp : param )
    {
        auto p = pp.second;
        if( p.name() == "Um" )
        {
            BOOST_CHECK_CLOSE( p.value(), 0.3, 10e-8);
            BOOST_CHECK_CLOSE( p.min(), 1e-4, 10e-8);
            BOOST_CHECK_EQUAL( p.max(), 10);
            BOOST_CHECK_EQUAL( p.description(), "Um desc" );
        }
        else if( p.name() == "H" )
        {
            BOOST_CHECK_CLOSE( p.value(), 0.41, 10e-8 );
            BOOST_CHECK_EQUAL( p.description(), "" );
        }
        else if( p.name() == "Te" )
        {
            BOOST_CHECK_CLOSE( p.value(), 293.1, 10e-8);
        }
#if 0
        Feel::cout << p.name() << std::endl
                   << "\tvalue : " << p.value() << std::endl
                   << "\tmin   : " << p.min() << std::endl
                   << "\tmax   : " << p.max() << std::endl
                   << "\tdesc  : " << p.description() << std::endl;
        if ( p.hasExpression() )
            Feel::cout << "\texpr  : " << p.expression() << std::endl;
#endif
    }
}

BOOST_AUTO_TEST_CASE( test_bc )
{
    std::set<std::string> markersRobinTest = {"mark1","mark2","mymarker","toto2","toto4","toto6","feelpp1_pA","feelpp1_pD","feelpp2_pA","feelpp2_pD"};

    ModelProperties model_props;
    model_props.setupFromFilenameOption();

    auto const& boundaryconditions = model_props.boundaryConditions();

    bool hasRobinBC = false;
    BOOST_CHECK( boundaryconditions.hasSection( "velocity" ) );
    for( auto const& [j_type,j_val] : boundaryconditions.section( "velocity" ).items() )
    {
        BOOST_CHECK( j_type == "Dirichlet" || j_type  == "Neumann" || j_type == "Robin" );
        BOOST_CHECK( j_val.is_object() );

        if ( j_type == "Robin" )
        {
            hasRobinBC = true;
            BOOST_CHECK( j_val.contains( "test" ) );
            auto const& j_robin_test = j_val.at( "test" );
            BOOST_CHECK( j_robin_test.contains( "markers" ) );
            ModelMarkers mmarkers;
            mmarkers.setup( j_robin_test.at( "markers" ) );

            for ( std::string const& _m : markersRobinTest )
                BOOST_CHECK( mmarkers.find( _m ) != mmarkers.end() );
        }
    }
    BOOST_CHECK( hasRobinBC );
}

BOOST_AUTO_TEST_CASE( test_bc2 )
{
    std::map<std::string,std::map<std::string,std::set<std::string>>> markers;
    markers["Dirichlet"] = { { "inflow", { "inflow" } }, { "wall", { "mark" } }, { "cylinder", { "mark1", "mark2" } } };
    markers["Neumann"] = { { "outlet", { "feelpp" } } };
    markers["Robin"] = { { "test", { "mark1","mark2","mymarker","toto2","toto4","toto6","feelpp1_pA","feelpp1_pD","feelpp2_pA","feelpp2_pD" } } };

    ModelProperties model_props;
    model_props.enableBoundaryConditions2();
    model_props.setupFromFilenameOption();

    auto boundaryconditions = model_props.boundaryConditions2();
    for( auto const& [bcTypeName,bcTypeData] : boundaryconditions.at("velocity") )
    {
        for( auto const& [bcName,cond] : bcTypeData )
        {
            auto const& checkMarkers = markers.at(bcTypeName).at(bcName);
            auto const& bcMarkers = cond.markers();
            for ( std::string const& _m : checkMarkers )
                BOOST_CHECK( bcMarkers.find( _m ) != bcMarkers.end() );

            if ( bcTypeName == "Robin" && bcName == "test" )
                BOOST_CHECK_EQUAL( cond.material(), "mycopper" );
        }
    }
}



BOOST_AUTO_TEST_CASE( test_postprocess )
{
    ModelProperties model_props;
    model_props.setupFromFilenameOption();

    std::map<std::string,std::tuple<double,std::set<std::string>> > ppStatExpected;
    for ( std::string const& i : std::vector<std::string>( {"A","B"} ) )
    {
        for ( int j : std::vector<int>( {3,5,7} ) )
        {
            std::get<0>( ppStatExpected[(boost::format("my_%1%_%2%_eval1")%i%j).str()] )=3.5*j;
            std::get<1>( ppStatExpected[(boost::format("my_%1%_%2%_eval1")%i%j).str()] )=
                { (boost::format("mat%1%%2%_x")%i%j).str(), (boost::format("mat%1%%2%_y")%i%j).str(), (boost::format("mat%1%%2%_z")%i%j).str() };
        }
        for ( auto [j,jstr] : std::vector<std::pair<int,std::string>>( { std::make_pair(3,"trois"), std::make_pair(5,"cinq"),std::make_pair(7,"sept") } ) )
        {
            std::get<0>( ppStatExpected[(boost::format("my_%1%_%2%_eval2")%i%jstr).str()] )=3.5*j;
            std::get<1>( ppStatExpected[(boost::format("my_%1%_%2%_eval2")%i%jstr).str()] )=
                { (boost::format("mat%1%%2%")%i%j).str() };
        }
    }

    std::map<std::string,std::tuple<double,std::set<std::string>> > ppStatRegistered;
    for (auto const& stat : model_props.postProcess().measuresStatistics("test") )
    {
        //std::cout << "stat.name() " << stat.name() << " " << stat.markers()  << std::endl;
        std::get<0>( ppStatRegistered[stat.name()] ) = stat.expr().exprScalar().evaluate()(0,0);
        std::get<1>( ppStatRegistered[stat.name()] ) = stat.markers();
    }

    for ( auto const& [name,statValues] : ppStatExpected )
    {
        auto itFindStat = ppStatRegistered.find( name );
        BOOST_CHECK( itFindStat != ppStatRegistered.end() );
        BOOST_CHECK_CLOSE( std::get<0>( itFindStat->second ), std::get<0>(statValues ), 1e-10 );
        for ( std::string const& marker : std::get<1>(statValues ) )
            BOOST_CHECK(  std::get<1>( itFindStat->second ).find( marker ) != std::get<1>( itFindStat->second ).end() );
    }

}

BOOST_AUTO_TEST_CASE( test_outputs )
{
    ModelProperties model_props;
    model_props.setupFromFilenameOption();

    auto outputs = model_props.outputs();
    for( auto const& out : outputs )
    {
        auto output = out.second;
        if( output.name() == "myoutput" )
        {
            BOOST_CHECK_EQUAL( output.type(), "average");
            BOOST_CHECK_EQUAL( output.markers().size(), 2);
            BOOST_CHECK_EQUAL( output.dim(), 3 );
        }
        else if( output.type() == "flux" )
        {
            BOOST_CHECK_EQUAL( output.type(), "flux");
            BOOST_CHECK_EQUAL( output.markers().size(), 1);
            BOOST_CHECK_EQUAL( output.dim(), 2 );
        }
#if 0
        std::cout << output;
#endif
    }
}

BOOST_AUTO_TEST_SUITE_END()
