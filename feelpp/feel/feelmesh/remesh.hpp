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
#pragma once

#include <mutex>
#include <mmg/libmmg.h>


namespace Feel {

/**
 * Class that handles remeshing in sequential and parallel 
 */
template<typename MeshType>
class Remesh
{
public:
    using mmg_mesh_t = MMG5_pMesh;
    using mesh_t = MeshType;
    using mesh_ptrtype = std::shared_ptr<MeshType>;
    using scalar_metric_t = typename Pch_type<mesh_t,1>::element_type;
    Remesh()
        :
        Remesh( nullptr )
    {}
    Remesh( std::shared_ptr<MeshType> const& mesh )
        :
        M_mesh( mesh ),
        M_mmg_mesh( nullptr ),
        M_mmg_sol( nullptr ),
        M_mmg_met( nullptr )
        {
            if constexpr ( dimension_v<MeshType> == 3 )
                MMG3D_Init_mesh(MMG5_ARG_start,
                                MMG5_ARG_ppMesh,&M_mmg_mesh,MMG5_ARG_ppMet,&M_mmg_sol,
                                MMG5_ARG_end);
            else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 3  )
                MMGS_Init_mesh(MMG5_ARG_start,
                               MMG5_ARG_ppMesh,&M_mmg_mesh,MMG5_ARG_ppMet,&M_mmg_sol,
                               MMG5_ARG_end);
            else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 2  )
                MMG2D_Init_mesh(MMG5_ARG_start,
                                MMG5_ARG_ppMesh,&M_mmg_mesh,MMG5_ARG_ppMet,&M_mmg_sol,
                                MMG5_ARG_end );
            
            this->mesh2Mmg();
            this->setParameters();
        }
    
    ~Remesh()
        {
            if constexpr ( dimension_v<MeshType> == 3 )
                MMG3D_Free_all(MMG5_ARG_start,
                               MMG5_ARG_ppMesh,&M_mmg_mesh,MMG5_ARG_ppMet,&M_mmg_sol,
                               MMG5_ARG_end);
            else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 3  )
                MMGS_Free_all(MMG5_ARG_start,
                              MMG5_ARG_ppMesh,&M_mmg_mesh,MMG5_ARG_ppMet,&M_mmg_sol,
                              MMG5_ARG_end);
            else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 2  )
                MMG2D_Free_all(MMG5_ARG_start,
                               MMG5_ARG_ppMesh,&M_mmg_mesh,MMG5_ARG_ppMet,&M_mmg_sol,
                               MMG5_ARG_end);
                
        }


    /**
     * set scalar metric
     */
    void setMetric( scalar_metric_t const& );

    /**
     * execute remesh task
     */
    mesh_ptrtype execute()
        {
            if constexpr ( dimension_v<MeshType> == 3 )
                MMG3D_mmg3dlib(M_mmg_mesh,M_mmg_sol);
            else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 3  )
                MMGS_mmgslib(M_mmg_mesh,M_mmg_sol);
            else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 2  )
                MMG2D_mmg2dlib(M_mmg_mesh,M_mmg_sol);
            auto r = this->mmg2Mesh();
            r->updateForUse();
            return r;
        }

    /**
     * transform a Feel++ \p mesh into an Mmg mesh
     */
    mmg_mesh_t mesh2Mmg( std::shared_ptr<MeshType> const& m_in );
    mmg_mesh_t mesh2Mmg() { return mesh2Mmg( M_mesh ); }
    
    /**
     * convert a Mmg mesh into a Feel++ \p mesh
     */
    mesh_ptrtype mmg2Mesh( mmg_mesh_t const& m_in );
    mesh_ptrtype mmg2Mesh() { return mmg2Mesh( M_mmg_mesh ); }

private:
    void setParameters();
    
private:

    std::shared_ptr<MeshType> M_mesh;
    MMG5_pMesh M_mmg_mesh;
    MMG5_pSol M_mmg_sol;
    MMG5_pSol M_mmg_met;

    std::unordered_map<int,int> pt_id;

    std::mutex mutex_;

};


template <typename MeshType>
Remesh<MeshType> remesher( std::shared_ptr<MeshType> const& m  )
{
    return Remesh<MeshType>{ m };
}

template<typename MeshType>
void Remesh<MeshType>::setMetric( scalar_metric_t const& m )
{
    if constexpr ( dimension_v<MeshType> == 3 )
    {
        if ( MMG3D_Set_solSize(M_mmg_mesh, M_mmg_sol,
                               MMG5_Vertex, M_mesh->numPoints(), MMG5_Scalar) != 1 )
        { 
            throw std::logic_error("Unable to allocate the metric array.");
        }
    }        
    else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 3  )
    {
        if ( MMGS_Set_solSize(M_mmg_mesh, M_mmg_sol,
                               MMG5_Vertex, M_mesh->numPoints(), MMG5_Scalar) != 1 )
        { 
            throw std::logic_error("Unable to allocate the metric array.");
        }
    }
    else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 2  )
    {
        if ( MMG2D_Set_solSize(M_mmg_mesh, M_mmg_sol,
                               MMG5_Vertex, M_mesh->numPoints(), MMG5_Scalar) != 1 )
        { 
            throw std::logic_error("Unable to allocate the metric array.");
        }
    }

    for( auto const& welt : M_mesh->elements() )
    {
        auto const& [key,elt] = boost::unwrap_ref( welt );
        for( auto const& ldof : m.functionSpace()->dof()->localDof( elt.id() ) )
        {
            size_type index = ldof.second.index();
            uint16_type local_dof = ldof.first.localDof();
            auto s = m( index );
            int pos = pt_id[elt.point( local_dof ).id()];
            if constexpr ( dimension_v<MeshType> == 3 )
            {
                if ( MMG3D_Set_scalarSol(M_mmg_sol,s,pos) != 1 )
                { 
                    throw std::logic_error("Unable to set metric");
                }
            }        
            else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 3  )
            {
                if ( MMGS_Set_scalarSol(M_mmg_sol,s,pos) != 1 )
                { 
                    throw std::logic_error("Unable to set metric");
                }
            }
            else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 2  )
            {
                if ( MMG2D_Set_scalarSol(M_mmg_sol,s,pos) != 1 )
                { 
                    throw std::logic_error("Unable to set metric");
                }
            }
        }
    }
}
template<typename MeshType>
MMG5_pMesh
Remesh<MeshType>::mesh2Mmg( std::shared_ptr<MeshType> const& m_in )
{
    int nVertices       = m_in->numPoints();
    int nTetrahedra     = ( dimension_v<MeshType> == 3 )?m_in->numElements():0;
    int nPrisms         = 0;
    int nTriangles      = ( dimension_v<MeshType> == 3 )?m_in->numFaces():m_in->numElements();
    int nQuadrilaterals = 0;
    int nEdges          = 0;
    
    if constexpr ( dimension_v<MeshType> == 3 )
    {
        if ( MMG3D_Set_meshSize(M_mmg_mesh,nVertices,nTetrahedra,nPrisms,nTriangles,
                                nQuadrilaterals,nEdges) != 1 )
        {
            throw std::logic_error( "Error in MMG3D_Set_meshSize" );
        }
    }
    else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 3  )
    {
        if ( MMGS_Set_meshSize(M_mmg_mesh,nVertices,nTriangles,nEdges) != 1 )
        {
            throw std::logic_error( "Error in MMGS_Set_meshSize" );
        }
    }
    else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 2 )
    {
        if ( MMG2D_Set_meshSize(M_mmg_mesh,nVertices,nTriangles,nEdges) != 1 )
        {
            throw std::logic_error( "Error in MMG2D_Set_meshSize" );
        }
    }

    
    pt_id.reserve( nVertices );
    int k = 1;
    for( auto const& [key,pt] : m_in->points() )
    {
        if constexpr ( dimension_v<MeshType> == 3 )
        {
            if ( MMG3D_Set_vertex(M_mmg_mesh, pt(0), pt(1), pt(2), pt.markerOr(0).value(), k) != 1 )
            {
                throw std::logic_error( "Error in MMG3D_Set_vertex" );
            }
            pt_id[pt.id()] = k++;
        }
        else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 3  )
        {
            if ( MMGS_Set_vertex(M_mmg_mesh, pt(0), pt(1), pt(2), pt.markerOr(0).value(), k) != 1 )
            {
                throw std::logic_error( "Error in MMGS_Set_vertex" );
            }
            pt_id[pt.id()] = k++;
        }
        else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 2  )
        {
            if ( MMG2D_Set_vertex(M_mmg_mesh, pt(0), pt(1), pt.markerOr(0).value(), k) != 1 )
            {
                throw std::logic_error( "Error in MMG2D_Set_vertex" );
            }
            pt_id[pt.id()] = k++;
        }
    }
    k = 1;
    for( auto const& welt : m_in->elements() )
    {
        auto const& [key,elt] = boost::unwrap_ref( welt );
        if constexpr ( dimension_v<MeshType> == 3 )
        {
            if ( MMG3D_Set_tetrahedron( M_mmg_mesh,
                                        pt_id[elt.point(0).id()],
                                        pt_id[elt.point(1).id()],
                                        pt_id[elt.point(2).id()],
                                        pt_id[elt.point(3).id()],
                                        elt.markerOr(0).value(), k ) != 1 )
            {
                throw std::logic_error( "Error in MMG3D_Set_tetrahedron" );
            }
        }
        else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 3  )
        {
            if ( MMGS_Set_triangle( M_mmg_mesh,
                                    pt_id[elt.point(0).id()],
                                    pt_id[elt.point(1).id()],
                                    pt_id[elt.point(2).id()],
                                    elt.markerOr(0).value(), k ) != 1 )
            {
                throw std::logic_error( "Error in MMGS_Set_triangle" );
            }
        }
        else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 2  )
        {
            if ( MMG2D_Set_triangle( M_mmg_mesh,
                                     pt_id[elt.point(0).id()],
                                     pt_id[elt.point(1).id()],
                                     pt_id[elt.point(2).id()],
                                     elt.markerOr(0).value(), k ) != 1 )
            {
                throw std::logic_error( "Error in MMG2D_Set_triangle" );
            }
        }
        k++;
    }

    if constexpr ( dimension_v<MeshType> == 3 )
    {
        int k = 1;
        for( auto const& wface : m_in->faces() )
        {
            auto const& [key,face] = boost::unwrap_ref( wface ); 
            if ( MMG3D_Set_triangle(M_mmg_mesh,
                                    pt_id[face.point(0).id()],
                                    pt_id[face.point(1).id()],
                                    pt_id[face.point(2).id()],
                                    face.markerOr(0).value(), k++ ) != 1 )
            {
            // TODO: hrow std::logic_error( "Error in MMG3D_Set_triangle" );
            }
        }
    }

    return M_mmg_mesh;
}

template<typename MeshType>
std::shared_ptr<MeshType> 
Remesh<MeshType>::mmg2Mesh( MMG5_pMesh const& mesh )
{
    int ier;

    int nVertices   = 0;
    int nTetrahedra = 0;
    int nTriangles  = 0;
    int nEdges      = 0;

    if ( MMG3D_Get_meshSize(M_mmg_mesh,&nVertices,&nTetrahedra,NULL,&nTriangles,NULL,
                            &nEdges) !=1 )
    { 
        ier = MMG5_STRONGFAILURE;
    }


    std::shared_ptr<MeshType> out = std::make_shared<MeshType>();
    int corner, required, tag;
    node_type n( mesh_t::nRealDim );
    for (int k = 1; k <= nVertices; k++)
    {
        if constexpr ( dimension_v<MeshType> == 3 )
        {
            if ( MMG3D_Get_vertex(mesh,&(n[0]),&(n[1]),&(n[2]),
                                  &(tag),&(corner),&(required)) != 1 )
            {
                cout << "Unable to get mesh vertex " << k << endl;
                ier = MMG5_STRONGFAILURE;
            }
        }
        else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 3 )
        {
            if ( MMGS_Get_vertex(mesh,&(n[0]),&(n[1]),&(n[2]),
                                 &(tag),&(corner),&(required)) != 1 )
            {
                cout << "Unable to get mesh vertex " << k << endl;
                ier = MMG5_STRONGFAILURE;
            }
        }
        if constexpr ( dimension_v<MeshType> == 2 )
        {
            
            if ( MMG2D_Get_vertex(mesh,&(n[0]),&(n[1]),
                                  &(tag),&(corner),&(required)) != 1 )
            {
                cout << "Unable to get mesh 2D vertex " << k << endl;
                ier = MMG5_STRONGFAILURE;
            }
        }
        using point_type = typename mesh_t::point_type;
        point_type pt( k, n );
        pt.setProcessIdInPartition( 0 );
        pt.setProcessId( 0 );
        pt.setMarker( tag );
        out->addPoint( pt );
    }

    if constexpr ( dimension_v<MeshType> == 3 )
    {
        for (int k=1; k<= nTetrahedra; k++ )
        {
            int iv[4], lab;
            if ( MMG3D_Get_tetrahedron(M_mmg_mesh,
                                       &(iv[0]),&(iv[1]),
                                       &(iv[2]),&(iv[3]),
                                       &(lab),&(required)) != 1 ) {
                std::cout << "Unable to get mesh tetra " << k << std::endl;
                ier = MMG5_STRONGFAILURE;
            }
            
            using element_type = typename mesh_t::element_type;
            element_type newElem;
            newElem.setMarker( lab );
            newElem.setProcessIdInPartition( 0 );
            newElem.setProcessId( 0 );
            for (int i=0; i<4; i++)
                newElem.setPoint( i, out->point( iv[i] ) );
            out->addElement( newElem,true );
        }
    }
    for ( int k=1; k<=nTriangles; k++ )
    {
        int iv[3], lab;
        if constexpr ( dimension_v<MeshType> == 3 )
        {
            if ( MMG3D_Get_triangle( M_mmg_mesh,
                                     &(iv[0]),&(iv[1]),&(iv[2]),
                                     &(lab),&(required)) != 1 )
            {
                std::cout << "Unable to get mesh triangle " << k << std::endl;
                ier = MMG5_STRONGFAILURE;
            }
            using face_type = typename mesh_t::face_type;
            face_type newElem;
            newElem.setMarker( lab );
            newElem.setProcessIdInPartition( 0 );
            newElem.setProcessId( 0 );
            for (int i=0; i<3; i++)
                newElem.setPoint( i, out->point( iv[i] ) );
            out->addFace( newElem );
        }
        if constexpr ( dimension_v<MeshType> == 2 &&  real_dimension_v<MeshType> == 2)
        {
            if ( MMG2D_Get_triangle( M_mmg_mesh,
                                     &(iv[0]),&(iv[1]),&(iv[2]),
                                     &(lab),&(required)) != 1 )
            {
                std::cout << "Unable to get mesh triangle " << k << std::endl;
                ier = MMG5_STRONGFAILURE;
            }
            using element_type = typename mesh_t::element_type;
            element_type newElem;
            newElem.setMarker( lab );
            newElem.setProcessIdInPartition( 0 );
            newElem.setProcessId( 0 );
            for (int i=0; i<3; i++)
                newElem.setPoint( i, out->point( iv[i] ) );
            out->addElement( newElem,true );
        }
        
    }

    std::cout << "vertices =" << nVertices << std::endl;
    std::cout << "tetrahedrons =" << nTetrahedra << std::endl;
    std::cout << "triangles =" << nTriangles << std::endl;
    std::cout << "Mesh" << out->numPoints() << " " << out->numElements() << " " << out->numFaces() << std::endl;

    return out;
}

template<typename MeshType>
void
Remesh<MeshType>::setParameters()
{
#if 0
    if constexpr ( dimension_v<MeshType> == 3 )
    {
        MMG3D_Set_iparameter(M_mmg_mesh,M_mmg_sol,MMG3D_IPARAM_verbose, value<MmgOption::Verbose>());
        MMG3D_Set_iparameter(M_mmg_mesh,M_mmg_sol,MMG3D_IPARAM_mem, value<MmgOption::Mem>());
        MMG3D_Set_dparameter(M_mmg_mesh,M_mmg_sol,MMG3D_DPARAM_hmin, value<MmgOption::Hmin>());
        MMG3D_Set_dparameter(M_mmg_mesh,M_mmg_sol,MMG3D_DPARAM_hmax, value<MmgOption::Hmax>());
        MMG3D_Set_dparameter(M_mmg_mesh,M_mmg_sol,MMG3D_DPARAM_hsiz, value<MmgOption::Hsiz>());
    }
    else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 3  )
    {
    }
    else if constexpr ( dimension_v<MeshType> == 2 && real_dimension_v<MeshType> == 2  )
    {

    }
#endif        

}
}


