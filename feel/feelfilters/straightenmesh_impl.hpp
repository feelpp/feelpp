/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2013-12-24

  Copyright (C) 2013 Feel++ Consortium

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
#if !defined(FEELPP_STRAIGHTENMESH_IMPL_HPP)
#define FEELPP_STRAIGHTENMESH_IMPL_HPP 1



#include <feel/feelcore/parameter.hpp>

#include <feel/feeldiscr/mesh.hpp>

#include <feel/feeldiscr/pchv.hpp>
#include <feel/feelmesh/meshmover.hpp>
#include <feel/feelfilters/detail/mesh.hpp>
#include <feel/feelvf/geometricdata.hpp>
#include <feel/feelvf/operators.hpp>
#include <feel/feelvf/integrate.hpp>
#include <feel/feelvf/projectors.hpp>
#include <feel/feelfilters/straightenmesh.hpp>

namespace Feel { namespace detail {
/// \cond DETAIL
template <typename ElementSpaceType>
void
straightenMeshUpdateEdgesOnBoundaryIsolated( ElementSpaceType & straightener, mpl::int_<0> /**/ )
{}
template <typename ElementSpaceType>
void
straightenMeshUpdateEdgesOnBoundaryIsolated( ElementSpaceType & straightener, mpl::int_<1> /**/ )
{}
template <typename ElementSpaceType>
void
straightenMeshUpdateEdgesOnBoundaryIsolated( ElementSpaceType & straightener, mpl::int_<2> /**/ )
{}
template <typename ElementSpaceType>
void
straightenMeshUpdateEdgesOnBoundaryIsolated( ElementSpaceType & straightener, mpl::int_<3> /**/ )
{
    typedef typename ElementSpaceType::functionspace_type space_type;
    typedef typename space_type::mesh_type mesh_type;
    typedef typename space_type::dof_type::fe_type fe_type;

    auto const ncdof = space_type::dof_type::nComponents;
    auto const dofshift = fe_type::nDofPerVertex*mesh_type::element_type::numVertices;

    auto mesh = straightener.mesh();
    auto const myrank = mesh->worldComm().localRank();

    std::set<size_type> edgeIdFoundToUpdate;

    // search if there are some edges on boundary which are the interface of two ghost boundary faces
    // in this case this edges must not be straigthen
    auto itedge = mesh->beginEdgeOnBoundary();
    auto const enedge = mesh->endEdgeOnBoundary();
    for ( ; itedge!=enedge ; ++itedge )
    {
        if ( itedge->numberOfProcGhost()==0 ) continue;

        auto const theedgeid = itedge->id();
        std::set<size_type> ghostFaceIdFoundOnBoundary;

        auto itprocghost=itedge->elementsGhost().begin();
        auto const enprocghost=itedge->elementsGhost().end();
        for ( ; itprocghost!=enprocghost ; ++itprocghost)
        {
            auto iteltghost = itprocghost->second.begin();
            auto const eneltghost = itprocghost->second.end();
            for ( ; iteltghost!=eneltghost ; ++iteltghost )
            {
                auto const& eltGhost = mesh->element(*iteltghost,itprocghost->first);
                for ( uint16_type f = 0 ; f < mesh_type::element_type::numTopologicalFaces ; ++f )
                {
                    auto const& theface = eltGhost.face(f);
                    if ( theface.isOnBoundary() )
                    {
                        bool findEdge=false;
                        for ( uint16_type e = 0; e < mesh_type::face_type::numEdges && !findEdge ; ++e )
                        {
                            if ( theface.edge(e).id() == theedgeid) { findEdge=true; ghostFaceIdFoundOnBoundary.insert(theface.id());}
                        }
                    }
                }
            }
        } // for ( ; itprocghost!=enprocghost ; ++itprocghost)

        // if 2 faces are find then the edge must not be straigten
        if (ghostFaceIdFoundOnBoundary.size()==2) edgeIdFoundToUpdate.insert(theedgeid);

    } // for ( ; itedge!=enedge ; ++itedge )

    // apply a null displacement on egdes founded (not straighten)
    if (edgeIdFoundToUpdate.size() > 0)
    {
        auto range = elements(mesh,EntityProcessType::ALL);
        auto iteltactif = range.template get<1>();
        auto const eneltactif = range.template get<2>();
        for ( ; iteltactif!=eneltactif ; ++iteltactif )
        {
            auto const& curElt = boost::unwrap_ref(*iteltactif);
            if ( !curElt.isOnBoundary() ) continue;

            for ( uint16_type e = 0; e < mesh_type::element_type::numEdges ; ++e )
            {
                if ( edgeIdFoundToUpdate.find(curElt.edge(e).id()) != edgeIdFoundToUpdate.end())
                {
                    // find edge
                    auto const idEltFind = curElt.id();
                    for ( uint16_type locdof = 0 ; locdof<fe_type::nDofPerEdge ; ++locdof )
                    {
                        auto const local_id = dofshift + e*fe_type::nDofPerEdge + locdof;
                        for ( uint16_type comp = 0; comp < ncdof; ++comp )
                        {
                            const size_type globdof = straightener.functionSpace()->dof()->localToGlobal( idEltFind, local_id, comp ).template get<0>();
                            straightener( globdof ) = 0;
                        }
                    }
                }
            }
        } // for ( ; iteltactif!=eneltactif ; ++iteltactif )
    } // if (edgeIdFoundToUpdate.size() > 0)


} // straightenMeshUpdateEdgesOnBoundaryIsolated

} // namespace detail
/// \endcond

/**
   \brief straighten the internal faces of a high order mesh

   \arg mesh mesh data structure
*/
template<typename MeshType>
boost::shared_ptr<MeshType>
straightenMesh( boost::shared_ptr<MeshType> mesh, WorldComm const& worldcomm, bool refine, bool save )
{
    typedef MeshType _mesh_type;
    typedef boost::shared_ptr<MeshType> _mesh_ptrtype;

    VLOG(1) << "straighten mesh of order " <<  _mesh_type::nOrder << " start";

    _mesh_ptrtype _mesh( mesh );

    using namespace vf;
    bool upExtendedElt = true;
    EntityProcessType entityProcess = (upExtendedElt)? EntityProcessType::ALL : EntityProcessType::LOCAL_ONLY;
    auto Xh = Pchv<_mesh_type::nOrder>( _mesh,upExtendedElt );
    auto xHo = vf::project( _space=Xh, _range=elements( mesh, entityProcess ), _expr=vf::P(), _geomap=GeomapStrategyType::GEOMAP_HO );
    auto xLo = vf::project( _space=Xh, _range=elements( mesh, entityProcess ), _expr=vf::P(), _geomap=GeomapStrategyType::GEOMAP_O1 );
    auto xHoBdy = vf::project( _space=Xh, _range=boundaryfaces( mesh, entityProcess ), _expr=vf::P(), _geomap=GeomapStrategyType::GEOMAP_HO );
    auto xLoBdy = vf::project( _space=Xh, _range=boundaryfaces( mesh, entityProcess ), _expr=vf::P(), _geomap=GeomapStrategyType::GEOMAP_O1 );

    auto straightener = Xh->element();
    straightener=( xLo-xHo )-( xLoBdy-xHoBdy );

    if (worldcomm.localSize()>1)
        Feel::detail::straightenMeshUpdateEdgesOnBoundaryIsolated( straightener,mpl::int_<_mesh_type::nDim>() );

    double norm_mean_value = integrate( _range=boundaryfaces( _mesh ), _expr=idv( straightener ) ).evaluate(true,worldcomm).norm();

    if ( norm_mean_value > 1e-12 )
        std::cout << "the straightening process may have failed\n"
                  << "norm of component-wise mean value of displacement on the boundary should be 0"
                  << "norm_mean_value: "  << norm_mean_value << "\n"
                  << "you should consider not using straightenMesh()\n"
                  << "\n";

#if 0
    boost::shared_ptr<Exporter<_mesh_type,_mesh_type::nOrder> > exporter;

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

    VLOG(1) << "straighten mesh of order " <<  _mesh_type::nOrder << " finish";

    return _mesh;
}



}
#endif /* FEELPP_STRAIGHTENMESH_IMPL_HPP */
