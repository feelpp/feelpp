/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

 This file is part of the Feel++ library

 Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org>
 Date: 2016-01-23

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
/**
   \file test_partitioner_metis.cpp
   \author Vincent Chabannes <vincent.chabannes@feelpp.org>
   \date 2016-01-23
 */

#define BOOST_TEST_MODULE test_partitioner_metis
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelpartition/partitionermetis.hpp>


FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( partitioner_metis )

BOOST_AUTO_TEST_CASE( partitioner_metis1 )
{
    using namespace Feel;
    std::vector<int> partIdsBuild;
    if ( Environment::isMasterRank() )
    {
        Metis::idx_t ncon = 1;

        Metis::idx_t options[METIS_NOPTIONS];
#if defined(FEELPP_USE_INTERNAL_METIS)
        Metis::Feel_METIS_SetDefaultOptions(options);
#else
        Metis::METIS_SetDefaultOptions(options);
#endif
        options[Metis::METIS_OPTION_NUMBERING]= 0;

        Metis::idx_t n = 4;
        Metis::idx_t nparts = 2;
        std::vector<Metis::idx_t> offset = { 0,2,4,6,8 };
        std::vector<Metis::idx_t> vals= {1,2,0,3,0,3,1,2};

        Metis::idx_t edgecut2;
        std::vector<Metis::idx_t> part(n);

#if defined(FEELPP_USE_INTERNAL_METIS)
        Metis::Feel_METIS_PartGraphRecursive(&n, &ncon, &offset[0], &vals[0], NULL, NULL,
                                             NULL, &nparts, NULL, NULL, options/*NULL*/,
                                             &edgecut2, &part[0]);
#else
        Metis::METIS_PartGraphRecursive(&n, &ncon, &offset[0], &vals[0], NULL, NULL,
                                             NULL, &nparts, NULL, NULL, options/*NULL*/,
                                             &edgecut2, &part[0]);
#endif

        std::set<int> partIdsBuildUnique( part.begin(), part.end() );
        partIdsBuild = std::vector<int>( partIdsBuildUnique.begin(), partIdsBuildUnique.end() );
    }
    mpi::broadcast( Environment::worldComm(), partIdsBuild, Environment::masterRank() );
    BOOST_CHECK( partIdsBuild.size() == 2 );
    BOOST_CHECK( *partIdsBuild.begin() == 0 );
    BOOST_CHECK( *partIdsBuild.rbegin() == 1 );
}

BOOST_AUTO_TEST_CASE( partitioner_metis2 )
{
    using namespace Feel;
    if ( Environment::isMasterRank() )
    {
        auto mesh = loadMesh(_mesh=new Mesh<Simplex<2>>( Environment::worldCommSeqPtr()),
                             _update=size_type(MESH_UPDATE_ELEMENTS_ADJACENCY|MESH_NO_UPDATE_MEASURES) );

        int numPartition = std::min((int)4, (int)mesh->numElements() );
        PartitionerMetis<decltype(mesh)> metis;
        metis.partition( mesh, numPartition );

        for (rank_type p=0;p<numPartition;++p)
        {
            auto rangeElements = mesh->elementsWithProcessId( p );
            int nElt = std::distance( std::get<0>( rangeElements ), std::get<1>( rangeElements ) );
            BOOST_CHECK( nElt > 0 );
        }
   }
}

BOOST_AUTO_TEST_SUITE_END()
