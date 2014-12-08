/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-07-15

  Copyright (C) 2014 Feel++ Consortium

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
#include <feel/feel.hpp>
#include <feel/feeldiscr/ned1h.hpp>

using namespace Feel;
struct DofEdgeInfo
{
    enum EdgeType
    { 
        EDGE_INTERIOR = 0, // edge in the interior
        EDGE_BOUNDARY, // edge on boundary 
        EDGE_BOUNDARY_VERTEX_1, // edge touches boundary with one vertex
        EDGE_BOUNDARY_VERTEX_2 // edge touches boundary with two vertices (but still is in the interior)
    };
    size_type index;
    size_type sign;
    EdgeType type;
    size_type dof_vertex_id1;
    size_type dof_vertex_id2;

};

int main(int argc, char**argv )
{
    
	po::options_description dofboundaryoptions( "dofboundary options" );
	dofboundaryoptions.add_options()
    ( "dof", po::value<int>()->default_value( 0 ), "global dof id" )
		;
	Environment env( _argc=argc, _argv=argv,
                     _desc= dofboundaryoptions,
                     _about=about(_name="dofboundary",
                                  _desc="print boundary dof information of a dof",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));

    auto mesh = loadMesh(_mesh=new Mesh<Simplex<3>>);
    auto Vh = Ned1h<0>( mesh );
    auto Xh = Pch<1>( mesh );

    std::vector<bool> done( Vh->nLocalDof(), false );
    // this vector store the ids
    std::vector<DofEdgeInfo> dof_edge_info( Vh->nLocalDof(), {invalid_uint16_type_value,1,DofEdgeInfo::EDGE_INTERIOR,invalid_size_type_value,invalid_size_type_value} );
    // loop over the multiset of global dof the dof table provides a relation
    // between local dof (element id + local dof id) and global dof. While the
    // global dof is unique, it is possibly associated to multiple local dof.
    // we can work here with the global view (global dof) of the relation or the
    // local view (local dof) of the relation.
    auto dofbegin = Vh->dof()->globalDof().first;
    auto dofend = Vh->dof()->globalDof().second;
    for( auto dofit = dofbegin; dofit != dofend; ++dofit)
    {
        auto const& dof = *dofit;
        if ( !done[ dof.first.index() ] )
        {
            // first print the local dof associated to the global dof 'dof'
            std::cout << "dof element id : " << dof.second.elementId() 
                      << " local dof : " << dof.second.localDof() << "\n";
            // note that local dof provides (for lowest order) the edge id for Nedelec
            
            auto const& edge = mesh->element(dof.second.elementId()).edge(dof.second.localDof());
            dof_edge_info[dof.first.index()].index = dof.first.index();
            dof_edge_info[dof.first.index()].sign = dof.first.sign();

            // get the vertex id of the end points
            auto const& pt1 = edge.point( 0 );
            auto const& pt2 = edge.point( 1 );
            size_type dofid1 = invalid_uint16_type_value;
            size_type dofid2 = invalid_uint16_type_value;
            // now find the vertex dof id in Xh associated to the point id from this edge
            // this is the id that is stored in the dof_edge info data structure
            for ( uint16_type i = 0; i < mesh->numLocalVertices(); ++i )
            {
                if ( mesh->element( dof.second.elementId() ).point( i ).id() == pt1.id() )
                {
                    dofid1 = Xh->dof()->localToGlobal( dof.second.elementId(), i, 0 ).index();
                }
                if ( mesh->element( dof.second.elementId() ).point( i ).id() == pt2.id() )
                {
                    dofid2 = Xh->dof()->localToGlobal( dof.second.elementId(), i, 0 ).index();
                }
                
            }

            if ( edge.isOnBoundary() )
            {
                dof_edge_info[dof.first.index()].type = DofEdgeInfo::EDGE_BOUNDARY;
                dof_edge_info[dof.first.index()].dof_vertex_id1 = dofid1;
                dof_edge_info[dof.first.index()].dof_vertex_id2 = dofid2;
            }
            if ( !edge.isOnBoundary() )
                {
                    
                    
                    //both points touch the boundary
                    if ( pt1.isOnBoundary() && pt2.isOnBoundary() )
                    {
                        dof_edge_info[dof.first.index()].type = DofEdgeInfo::EDGE_BOUNDARY_VERTEX_2;   
                        dof_edge_info[dof.first.index()].dof_vertex_id1 = dofid1;
                        dof_edge_info[dof.first.index()].dof_vertex_id2 = dofid2;
                        CHECK( dofid1 != invalid_size_type_value ) << "Invalid dof vertex id1";
                        CHECK( dofid2 != invalid_size_type_value ) << "Invalid dof vertex id2";
                    }
                    // one of the end points touch the boundary
                    else if ( pt1.isOnBoundary()  )
                    {
                        dof_edge_info[dof.first.index()].type = DofEdgeInfo::EDGE_BOUNDARY_VERTEX_1;   
                        dof_edge_info[dof.first.index()].dof_vertex_id1 = dofid1;
                        dof_edge_info[dof.first.index()].dof_vertex_id2 = invalid_size_type_value;
                        CHECK( dofid1 != invalid_size_type_value ) << "Invalid dof vertex id1";
                    }
                    else if ( pt2.isOnBoundary()  )
                    {
                        dof_edge_info[dof.first.index()].type = DofEdgeInfo::EDGE_BOUNDARY_VERTEX_1;   
                        dof_edge_info[dof.first.index()].dof_vertex_id1 = dofid2;
                        dof_edge_info[dof.first.index()].dof_vertex_id2 = invalid_size_type_value;
                        CHECK( dofid2 != invalid_size_type_value ) << "Invalid dof vertex id1";
                    }
                    else
                    {
                        dof_edge_info[dof.first.index()].type = DofEdgeInfo::EDGE_INTERIOR;   
                        dof_edge_info[dof.first.index()].dof_vertex_id1 = invalid_size_type_value;
                        dof_edge_info[dof.first.index()].dof_vertex_id2 = invalid_size_type_value;
                    }   
                    
                }
            std::cout << "  - edge id " << dof_edge_info[dof.first.index()].index << " type " << dof_edge_info[dof.first.index()].type << " sign " << dof_edge_info[dof.first.index()].sign << "\n";
        }
        done[dof.first.index()] = true;
            
    }
    return 0;

}
