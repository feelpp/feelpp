/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2013-12-24

  Copyright (C) 2013-2016 Feel++ Consortium

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
   \file straightenmesh.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2013-12-24
 */
#if !defined( FEELPP_STRAIGHTENMESH_IMPL_HPP )
#define FEELPP_STRAIGHTENMESH_IMPL_HPP 1

#include <feel/feelcore/parameter.hpp>

#include <feel/feeldiscr/mesh.hpp>

#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelfilters/detail/mesh.hpp>
#include <feel/feelfilters/straightenmesh.hpp>
#include <feel/feelmesh/meshmover.hpp>
#include <feel/feelvf/geometricdata.hpp>
#include <feel/feelvf/integrate.hpp>
#include <feel/feelvf/operators.hpp>
#include <feel/feelvf/projectors.hpp>

namespace Feel
{
/**
   \brief straighten the internal faces of a high order mesh

   \arg mesh mesh data structure
*/
template <typename MeshType>
std::shared_ptr<MeshType>
straightenMesh( std::shared_ptr<MeshType> mesh, worldcomm_ptr_t const& worldcomm, bool refine, bool save )
{
    if constexpr ( MeshType::nOrder > 1 )
    {
        typedef MeshType _mesh_type;
        typedef std::shared_ptr<MeshType> _mesh_ptrtype;
        using size_type = typename MeshType::size_type;
        VLOG( 1 ) << "straighten mesh of order " << _mesh_type::nOrder << " start";

        _mesh_ptrtype _mesh( mesh );

        using namespace vf;
        bool upExtendedElt = true;
        EntityProcessType entityProcess = ( upExtendedElt ) ? EntityProcessType::ALL : EntityProcessType::LOCAL_ONLY;
        auto Xh = Pchv<_mesh_type::nOrder>( _mesh, upExtendedElt );
        auto xHo = vf::project( _space = Xh, _range = elements( mesh, entityProcess ), _expr = vf::P(), _geomap = GeomapStrategyType::GEOMAP_HO );
        auto xLo = vf::project( _space = Xh, _range = elements( mesh, entityProcess ), _expr = vf::P(), _geomap = GeomapStrategyType::GEOMAP_O1 );
        auto xHoBdy = vf::project( _space = Xh, _range = boundaryfaces( mesh, entityProcess ), _expr = vf::P(), _geomap = GeomapStrategyType::GEOMAP_HO );
        auto xLoBdy = vf::project( _space = Xh, _range = boundaryfaces( mesh, entityProcess ), _expr = vf::P(), _geomap = GeomapStrategyType::GEOMAP_O1 );

        std::set<size_type> dofsWithValueImposed;
        for ( auto const& faceWrap : boundaryfaces( mesh, entityProcess ) )
        {
            auto const& face = unwrap_ref( faceWrap );
            auto facedof = Xh->dof()->faceLocalDof( face.id() );
            for ( auto it = facedof.first, en = facedof.second; it != en; ++it )
                dofsWithValueImposed.insert( it->index() );
        }
        sync( xHoBdy, "=", dofsWithValueImposed );
        sync( xLoBdy, "=", dofsWithValueImposed );

        auto straightener = Xh->element();
        straightener = ( xLo - xHo ) - ( xLoBdy - xHoBdy );

#if 0
    if (worldcomm->localSize()>1)
        Feel::detail::straightenMeshUpdateEdgesOnBoundaryIsolated( straightener,mpl::int_<_mesh_type::nDim>() );
#endif
        double norm_mean_value = integrate( _range = boundaryfaces( _mesh ), _expr = idv( straightener ) ).evaluate( parallelEvaluation ).norm();

        if ( norm_mean_value > 1e-12 )
            std::cout << "the straightening process may have failed\n"
                      << "norm of component-wise mean value of displacement on the boundary should be 0"
                      << "norm_mean_value: " << norm_mean_value << "\n"
                      << "you should consider not using straightenMesh()\n"
                      << "\n";

#if 0
    std::shared_ptr<Exporter<_mesh_type,_mesh_type::nOrder> > exporter;

    if ( save )
    {
        exporter = Exporter<_mesh_type,_mesh_type::nOrder>::New( "gmsh"/*test_app->vm()*/, "straightener" );
        exporter->step( 0 )->setMesh( _mesh );
        exporter->step( 0 )->add( "xHo", xHo );
        exporter->step( 0 )->add( "xLo", xLo );
        exporter->step( 0 )->add( "xHoBdy", xHoBdy );
        exporter->step( 0 )->add( "xLoBdy", xLoBdy );
        exporter->step( 0 )->add( "straightener", straightener );
        exporter->save();
    }
#endif

        MeshMover<_mesh_type> meshmove;
        meshmove.apply( _mesh, straightener );

        VLOG( 1 ) << "straighten mesh of order " << _mesh_type::nOrder << " finish";

        return _mesh;
    }
    else
        return mesh;
}

} // namespace Feel
#endif /* FEELPP_STRAIGHTENMESH_IMPL_HPP */
