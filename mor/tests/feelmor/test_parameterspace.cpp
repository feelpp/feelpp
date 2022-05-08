/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date: 2016-08-07

 Copyright (C) 2014-2016 Feel++ Consortium

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

#define BOOST_TEST_MODULE parameterspace testsuite
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelmor/parameterspace.hpp>

using namespace Feel;

inline
po::options_description
makeOptions()
{
    po::options_description modelopt( "test parameter space options" );
    modelopt.add_options()
        ("json_filename" , Feel::po::value<std::string>()->default_value( "$cfgdir/parameterspace.json" ), "json file" )
        ;
    return  modelopt.add( Feel::feel_options() ) ;
}

inline
AboutData
makeAbout()
{
    AboutData about( "test_parameterspace",
                     "test_parameterspace" );
    return about;

}

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() )

BOOST_AUTO_TEST_SUITE( parameterspace )

BOOST_AUTO_TEST_CASE( test1 )
{
    using namespace Feel;
    typedef ParameterSpace</*4*/> parameterspace_type;
    typedef typename parameterspace_type::element_type parameter_type;

    auto muspace = parameterspace_type::New(8);
    muspace->setDimension(4);

    auto muMin = muspace->element();
    muMin << 2,3,-4,5;
    parameter_type muMax(muspace);
    muMax(0)=6;muMax(1)=7;muMax(2)=8;muMax(3)=9;
    muspace->setMin( muMin );
    muspace->setMax( muMax );
    muspace->setParameterName( 0, "myparamA" );
    muspace->setParameterName( 1, "myparamB" );
    muspace->setParameterName( 2, "myparamC" );
    muspace->setParameterName( 3, "myparamD" );

    auto mysampling = muspace->sampling();
    // equidistribute sampling (identical for all process)
    mysampling->equidistribute( 120, true );
    int nSamples = mysampling->size();
    BOOST_CHECK( nSamples <= 120 );
    for ( int k=0;k<nSamples;++k )
        for ( uint16_type d=0;d<muspace->dimension();++d)
            BOOST_CHECK( (*mysampling)[k](d) >= muMin(d) && (*mysampling)[k](d) <= muMax(d) );
    // equidistribute sampling (distributed on world)
    mysampling->equidistribute( 120, false );
    BOOST_CHECK( mysampling->size() <= nSamples );
    for ( int k=0;k<mysampling->size();++k )
        for ( uint16_type d=0;d<muspace->dimension();++d)
            BOOST_CHECK( (*mysampling)[k](d) >= muMin(d) && (*mysampling)[k](d) <= muMax(d) );
    int nSamples2 = 0;
    mpi::all_reduce( Environment::worldComm(), int(mysampling->size()), nSamples2, std::plus<int>() );
    BOOST_CHECK( nSamples == nSamples2 );

    // equidistributeProduct sampling (identical for all process)
    std::vector<size_type> nSampleDir = { 12,5,8,3 };
    mysampling->equidistributeProduct( nSampleDir, true );
    nSamples = mysampling->size();
    BOOST_CHECK( nSamples == 12*5*8*3 );
    for ( int k=0;k<nSamples;++k )
        for ( uint16_type d=0;d<muspace->dimension();++d)
            BOOST_CHECK( (*mysampling)[k](d) >= muMin(d) && (*mysampling)[k](d) <= muMax(d) );
    // equidistributeProduct sampling (distributed on world)
    mysampling->equidistributeProduct( nSampleDir, false );
    BOOST_CHECK( mysampling->size() <= nSamples );
    for ( int k=0;k<mysampling->size();++k )
        for ( uint16_type d=0;d<muspace->dimension();++d)
            BOOST_CHECK( (*mysampling)[k](d) >= muMin(d) && (*mysampling)[k](d) <= muMax(d) );
    mpi::all_reduce( Environment::worldComm(), int(mysampling->size()), nSamples2, std::plus<int>() );
    BOOST_CHECK( nSamples == nSamples2 );

    // logEquidistribute sampling (identical for all process)
    mysampling->logEquidistribute( 230 );
    nSamples = mysampling->size();
    BOOST_CHECK( nSamples <= 230 );
    for ( int k=0;k<nSamples;++k )
        for ( uint16_type d=0;d<muspace->dimension();++d)
            BOOST_CHECK( (*mysampling)[k](d) >= muMin(d) && (*mysampling)[k](d) <= muMax(d) );
    // logEquidistribute sampling (distributed on world)
    mysampling->equidistribute( 230, false );
    BOOST_CHECK( mysampling->size() <= nSamples );
    for ( int k=0;k<mysampling->size();++k )
        for ( uint16_type d=0;d<muspace->dimension();++d)
            BOOST_CHECK( (*mysampling)[k](d) >= muMin(d) && (*mysampling)[k](d) <= muMax(d) );
    mpi::all_reduce( Environment::worldComm(), int(mysampling->size()), nSamples2, std::plus<int>() );
    BOOST_CHECK( nSamples == nSamples2 );

    // logEquidistributeProduct sampling (identical for all process)
    mysampling->logEquidistributeProduct( nSampleDir, true );
    nSamples = mysampling->size();
    BOOST_CHECK( nSamples == 12*5*8*3 );
    for ( int k=0;k<nSamples;++k )
        for ( uint16_type d=0;d<muspace->dimension();++d)
            BOOST_CHECK( (*mysampling)[k](d) >= muMin(d) && (*mysampling)[k](d) <= muMax(d) );
    // equidistributeProduct sampling (distributed on world)
    mysampling->logEquidistributeProduct( nSampleDir, false );
    BOOST_CHECK( mysampling->size() <= nSamples );
    for ( int k=0;k<mysampling->size();++k )
        for ( uint16_type d=0;d<muspace->dimension();++d)
            BOOST_CHECK( (*mysampling)[k](d) >= muMin(d) && (*mysampling)[k](d) <= muMax(d) );
    mpi::all_reduce( Environment::worldComm(), int(mysampling->size()), nSamples2, std::plus<int>() );
    BOOST_CHECK( nSamples == nSamples2 );

    // randomize sampling (identical for all process)
    mysampling->randomize( 168, true, "", false );
    nSamples = mysampling->size();
    BOOST_CHECK( nSamples == 168 );
    for ( int k=0;k<nSamples;++k )
        for ( uint16_type d=0;d<muspace->dimension();++d)
            BOOST_CHECK( (*mysampling)[k](d) >= muMin(d) && (*mysampling)[k](d) <= muMax(d) );
    // log-randomize sampling (distributed on world)
    mysampling->randomize( 168, false, "", true );
    BOOST_CHECK( mysampling->size() <= nSamples );
    for ( int k=0;k<mysampling->size();++k )
        for ( uint16_type d=0;d<muspace->dimension();++d)
            BOOST_CHECK( (*mysampling)[k](d) >= muMin(d) && (*mysampling)[k](d) <= muMax(d) );
    mpi::all_reduce( Environment::worldComm(), int(mysampling->size()), nSamples2, std::plus<int>() );
    BOOST_CHECK( nSamples == nSamples2 );

    // write sampling in json file
    mysampling->saveJson("mysampling.json");

    // create another sampling
    mysampling->equidistribute( 1249 );
    auto mysampling2 = muspace->sampling();
    mysampling2->setSuperSampling( mysampling );
    for ( int k=0;k<mysampling->size();++k )
    {
        if ( k%3 != 0 ) continue;
        mysampling2->push_back( (*mysampling)[k] );
    }

    // test complement
    auto mysamplingComplement = mysampling2->complement();
    BOOST_CHECK( mysampling->size() == ( mysampling2->size() + mysamplingComplement->size() ) );

#if defined(FEELPP_HAS_ANN_H)
    // test searchNearestNeighbors
    auto muSearch = muspace->element();
    for ( uint16_type d=0;d<muspace->dimension();++d)
        BOOST_CHECK( muSearch(d) >= muMin(d) && muSearch(d) <= muMax(d) );
    std::vector<int> indexNeighbors;
    auto mysamplingNeighbors = mysampling->searchNearestNeighbors( muSearch,5,indexNeighbors );
    BOOST_CHECK( mysamplingNeighbors->size() == 5 );
    mysamplingNeighbors->saveJson("mysamplingNeighbors.json");
#endif
    // min/max (identical for all process)
    std::vector<size_type> nSampleDir2 = { 5, 5, 13, 5 };
    mysampling->equidistributeProduct( nSampleDir2, true );
    auto minRes = mysampling->min( true );
    auto maxRes = mysampling->max( true );
    auto const& muMinRes = boost::get<0>( minRes );
    auto const& muMaxRes = boost::get<0>( maxRes );
    int indexMinRes = boost::get<1>( minRes );
    int indexMaxRes = boost::get<1>( maxRes );
    BOOST_CHECK( muMinRes(0) == 2 && muMinRes(1) == 3 && muMinRes(2) == 0 && muMinRes(3) == 5 );
    BOOST_CHECK( muMaxRes(0) == 6 && muMaxRes(1) == 7 && muMaxRes(2) == 8 && muMaxRes(3) == 9 );
    BOOST_CHECK( muMinRes == mysampling->at( indexMinRes ) );
    BOOST_CHECK( muMaxRes == mysampling->at( indexMaxRes ) );
    // min/max (distributed on world)
    mysampling->equidistributeProduct( nSampleDir2, false );
    auto minRes2 = mysampling->min( true );
    auto maxRes2 = mysampling->max( true );
    auto const& muMinRes2 = boost::get<0>( minRes2 );
    auto const& muMaxRes2 = boost::get<0>( maxRes2 );
    BOOST_CHECK( muMinRes2(0) == 2 && muMinRes2(1) == 3 && muMinRes2(2) == 0 && muMinRes2(3) == 5 );
    BOOST_CHECK( muMaxRes2(0) == 6 && muMaxRes2(1) == 7 && muMaxRes2(2) == 8 && muMaxRes2(3) == 9 );
}
BOOST_AUTO_TEST_CASE( test2 )
{
    using namespace Feel;
    typedef ParameterSpace</*4*/> parameterspace_type;
    typedef typename parameterspace_type::element_type parameter_type;

    auto muspace = parameterspace_type::New(4);
    auto muMin = muspace->element();
    muMin << 2,3,-4,3;
    auto muMax = muspace->element();
    muMax << 8,7,-2,5;
    muspace->setMin( muMin );
    muspace->setMax( muMax );
    muspace->setParameterName( 0, "myparamA" );
    muspace->setParameterName( 1, "myparamB" );
    muspace->setParameterName( 2, "myparamC" );
    muspace->setParameterName( 3, "myparamD" );
    muspace->saveJson("test2_muspace.json");

    auto muspaceReloaded = parameterspace_type::New("test2_muspace.json");
    muspaceReloaded->saveJson("test2_muspaceReloaded.json");
    BOOST_CHECK( muspace->dimension() == muspaceReloaded->dimension() );
    for (int d=0;d<muspace->dimension();++d )
    {
        BOOST_CHECK_CLOSE( muspace->min()(d),muspaceReloaded->min()(d),1e-9 );
        BOOST_CHECK( muspace->parameterName(d) == muspaceReloaded->parameterName(d) );
    }
}
BOOST_AUTO_TEST_CASE( test3 )
{
    using namespace Feel;
    using parameterspace_type = ParameterSpace<>;
    namespace pt =  boost::property_tree;

    auto json = removeComments(readFromFile(Environment::expand(soption("json_filename"))));
    std::istringstream istr( json );
    nl::json p;
    istr >> p;
    //pt::ptree p;
    //pt::read_json(istr, p);
    auto parameters = CRBModelParameters();
    parameters.setPTree(p);
    auto Dmu = parameterspace_type::New(parameters);
    BOOST_CHECK( Dmu->dimension() == 2 );
    // auto mubar = Dmu->mubar();
    // BOOST_CHECK( mubar.parameterNamed("p1") == 2 );
    // BOOST_CHECK( mubar.parameterNamed("p2") == -1 );
    auto mumin = Dmu->min();
    BOOST_CHECK( mumin.parameterNamed("p1") == 1 );
    BOOST_CHECK( mumin.parameterNamed("p2") == -10 );
    auto mumax = Dmu->max();
    BOOST_CHECK( mumax.parameterNamed("p1") == 3 );
    BOOST_CHECK( mumax.parameterNamed("p2") == 10 );
}
BOOST_AUTO_TEST_SUITE_END()
