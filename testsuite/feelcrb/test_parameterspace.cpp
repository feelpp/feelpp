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
#include <testsuite/testsuite.hpp>

#include <feel/feelcrb/parameterspace.hpp>

FEELPP_ENVIRONMENT_NO_OPTIONS

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
    // equidistribute sampling
    int nSampling = 120;
    mysampling->equidistribute( nSampling );
    BOOST_CHECK( mysampling->size() == nSampling );
    for ( int k=0;k<nSampling;++k )
        for ( uint16_type d=0;d<muspace->dimension();++d)
            BOOST_CHECK( (*mysampling)[k](d) >= muMin(d) && (*mysampling)[k](d) <= muMax(d) );

    // logEquidistribute sampling
    nSampling = 230;
    mysampling->logEquidistribute( nSampling );
    BOOST_CHECK( mysampling->size() == nSampling );
    for ( int k=0;k<nSampling;++k )
        for ( uint16_type d=0;d<muspace->dimension();++d)
            BOOST_CHECK( (*mysampling)[k](d) >= muMin(d) && (*mysampling)[k](d) <= muMax(d) );

    // randomize sampling
    nSampling = 168;
    mysampling->randomize( nSampling );
    BOOST_CHECK( mysampling->size() == nSampling );
    for ( int k=0;k<nSampling;++k )
        for ( uint16_type d=0;d<muspace->dimension();++d)
            BOOST_CHECK( (*mysampling)[k](d) >= muMin(d) && (*mysampling)[k](d) <= muMax(d) );

    // write sampling in json file
    mysampling->saveJson("mysampling.json");

    // create another sampling
    nSampling = 1249;
    mysampling->equidistribute( nSampling );
    auto mysampling2 = muspace->sampling();
    mysampling2->setSuperSampling( mysampling );
    for ( int k=0;k<nSampling;++k )
    {
        if ( k%3 != 0 ) continue;
        mysampling2->push_back( (*mysampling)[k] );
    }

    // test complement
    auto mysamplingComplement = mysampling2->complement();
    BOOST_CHECK( nSampling == ( mysampling2->size() + mysamplingComplement->size() ) );

    // test searchNearestNeighbors
    auto muSearch = muspace->element();
    for ( uint16_type d=0;d<muspace->dimension();++d)
        BOOST_CHECK( muSearch(d) >= muMin(d) && muSearch(d) <= muMax(d) );
    std::vector<int> indexNeighbors;
    auto mysamplingNeighbors = mysampling->searchNearestNeighbors( muSearch,5,indexNeighbors );
    BOOST_CHECK( mysamplingNeighbors->size() == 5 );
    mysamplingNeighbors->saveJson("mysamplingNeighbors.json");

}
BOOST_AUTO_TEST_SUITE_END()
