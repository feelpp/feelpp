/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Gilles Steiner <gilles.steiner@epfl.ch>
             Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2004-11-18

  Copyright (C) 2004 EPFL

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file filterfromvtk.hpp
   \author Gilles Steiner <gilles.steiner@epfl.ch>
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2004-11-18
 */
#ifndef __filter_H
#define __filter_H 1

#include <feel/feelcore/visitor.hpp>
//#include <feel/feelmesh/mesh1d.hpp>
#include <feel/feelmesh/mesh2d.hpp>
#include <feel/feelmesh/mesh3d.hpp>
#include <feel/feeldiscr/mesh.hpp>

#include <feel/feelcore/feel.hpp>


#if defined(FEELPP_HAS_VTK)
// Vtk header files
//#include "vtkPointSet.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkCell.h"

#endif /* FEELPP_HAS_VTK */

namespace Feel
{

/**
 * \class FilterFromVtk
 *
 * Converts Mesh data structure from Vtk library to Feel Mesh type.
 *
 * \author Gilles Steiner <gilles.steiner@epfl.ch>
 * \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 */
template<typename MeshType>
class FilterFromVtk
    :
public VisitorBase,
public Visitor<MeshType>
{

public:

    static const uint16_type nDim = MeshType::nDim;
    BOOST_STATIC_ASSERT( nDim == 2 || nDim == 3 );

    /** @name Typedefs
     */
    //@{
    typedef MeshType mesh_type;
    typedef typename mesh_type::point_type point_type;
    typedef typename point_type::node_type node_type;
    typedef typename mesh_type::edge_type edge_type;
    typedef typename mesh_type::face_type face_type;
    typedef typename mesh_type::element_type element_type;
#if defined(FEELPP_HAS_VTK)

    typedef typename mpl::if_<mpl::equal_to<mpl::int_<MeshType::nDim>,mpl::int_<2> >,
            mpl::identity<vtkPolyData>,
            mpl::identity<vtkUnstructuredGrid> >::type::type vtkmesh_type;


#endif /* FEELPP_HAS_VTK */
    //@}

    /** @name Constructors, destructor
     */
    //@{

#if defined(FEELPP_HAS_VTK)
    FilterFromVtk( vtkmesh_type* __vtkmesh )
    {
        _M_vtkmesh = vtkmesh_type::New();
        _M_vtkmesh->CopyStructure( __vtkmesh );
    }
#else
    FilterFromVtk()
    {}
#endif /* FEELPP_HAS_VTK */

    ~FilterFromVtk()
    {
#if defined(FEELPP_HAS_VTK)
        _M_vtkmesh->Delete();
#endif /* FEELPP_HAS_VTK */
    }

#if defined(FEELPP_HAS_VTK)
    vtkmesh_type* getVtkMesh()
    {
        return _M_vtkmesh;
    }
#endif /* FEELPP_HAS_VTK */

    //@}

    /** @name  Methods
     */
    //@{

    void visit( mesh_type* mesh )
    {
        visit( mesh, mpl::int_<nDim>() );
    }


    //@}

protected :

#if defined(FEELPP_HAS_VTK)
    vtkmesh_type* _M_vtkmesh;
#endif /* FEELPP_HAS_VTK */

private:

    void visit( mesh_type* mesh, mpl::int_<2> );
    //void visit( mesh_type* mesh, mpl::int_<3> );
};

/**
 * \class FilterFromVtk3D
 *
 * Converts Mesh data structure from Vtk library to Feel Mesh type.
 *
 * \author Gilles Steiner <gilles.steiner@epfl.ch>
 * \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 */
template<typename MeshType>
class FilterFromVtk3D
    :
public VisitorBase,
public Visitor<MeshType>
{

public:

    static const uint16_type nDim = MeshType::nDim;
    BOOST_STATIC_ASSERT( nDim == 2 || nDim == 3 );

    /** @name Typedefs
     */
    //@{
    typedef MeshType mesh_type;
    typedef typename mesh_type::point_type point_type;
    typedef typename point_type::node_type node_type;
    typedef typename mesh_type::edge_type edge_type;
    typedef typename mesh_type::face_type face_type;
    typedef typename mesh_type::element_type element_type;
#if defined(FEELPP_HAS_VTK)
    typedef typename mpl::if_<mpl::equal_to<mpl::int_<MeshType::nDim>,mpl::int_<2> >,
            mpl::identity<vtkPolyData>,
            mpl::identity<vtkUnstructuredGrid> >::type::type vtkmesh_type;

#endif /* FEELPP_HAS_VTK */
    //@}

    /** @name Constructors, destructor
     */
    //@{

#if defined(FEELPP_HAS_VTK)
    FilterFromVtk3D( vtkmesh_type* __vtkmesh )
    {
        _M_vtkmesh = vtkmesh_type::New();
        _M_vtkmesh->CopyStructure( __vtkmesh );
    }
#else
    FilterFromVtk3D()
    {}
#endif /* FEELPP_HAS_VTK */

    ~FilterFromVtk3D()
    {
#if defined(FEELPP_HAS_VTK)
        _M_vtkmesh->Delete();
#endif /* FEELPP_HAS_VTK */
    }

#if defined(FEELPP_HAS_VTK)
    vtkmesh_type* getVtkMesh()
    {
        return _M_vtkmesh;
    }
#endif /* FEELPP_HAS_VTK */

    //@}

    /** @name  Methods
     */
    //@{

    void visit( mesh_type* mesh )
    {
        visit( mesh, mpl::int_<nDim>() );
    }


    //@}

protected :

#if defined(FEELPP_HAS_VTK)
    vtkmesh_type* _M_vtkmesh;
#endif /* FEELPP_HAS_VTK */

private:

    //void visit( mesh_type* mesh, mpl::int_<2> );
    void visit( mesh_type* mesh, mpl::int_<3> );
};

template<typename MeshType>
void
FilterFromVtk<MeshType>::visit( mesh_type* mesh, mpl::int_<2> )
{
    detail::ignore_unused_variable_warning( mesh );
#if defined(FEELPP_HAS_VTK)
    //  std::cout <<"Start of mesh conversion !" << std::endl;

    vtkPolyData * _vtkMesh = this->getVtkMesh();

    uint __n = _vtkMesh->GetNumberOfPoints(); // Number of nodes

    Debug( 8099 ) <<"Number of points : "<< __n << "\n";

    uint __nele = _vtkMesh->GetNumberOfPolys(); // Number of elements

    Debug( 8099 ) <<"Number of elements : "<< __nele << "\n";

    // add the points to the mesh

    for ( uint __i = 0; __i < __n; ++__i )
    {
        node_type __nd( 2 );
        __nd[0] = _vtkMesh->GetPoint( __i )[0];
        __nd[1] = _vtkMesh->GetPoint( __i )[1];
        point_type __pt( __i,__nd, false );
        __pt.marker() = 0;

        if ( __nd[0] == -1 || __nd[1] == -1 || __nd[0] + __nd[1] == 0 )

        {
            __pt.setOnBoundary( true );

        }

        else
        {
            __pt.setOnBoundary( false );
        }

        mesh->addPoint( __pt );
    }


    size_type n_faces = 0;

    // Add Boundary faces

    face_type* pf0 = new face_type;

    pf0->setMarker( 1 );
    pf0->setPoint( 0, mesh->point( 1 ) );
    pf0->setPoint( 1, mesh->point( 2 ) );

    pf0->setId( n_faces++ );
    pf0->setOnBoundary( true );
    mesh->addFace( *pf0 );

    delete pf0;

    face_type* pf1 = new face_type;

    pf1->setMarker( 1 );
    pf1->setPoint( 0, mesh->point( 2 ) );
    pf1->setPoint( 1, mesh->point( 0 ) );

    pf1->setId( n_faces++ );
    pf1->setOnBoundary( true );
    mesh->addFace( *pf1 );

    delete pf1;

    face_type* pf2 = new face_type;

    pf2->setMarker( 1 );
    pf2->setPoint( 0, mesh->point( 0 ) );
    pf2->setPoint( 1, mesh->point( 1 ) );

    pf2->setId( n_faces++ );
    pf2->setOnBoundary( true );
    mesh->addFace( *pf2 );

    delete pf2;

    // add the elements to the mesh

    for ( uint __i = 0; __i < __nele; ++__i )
    {
        Debug( 8099 ) << "[FilterFromVtk] element " << __i << "\n";
        // Here we only have triangular elements of order 1

        element_type * pf = new element_type;

        pf->setId( __i );
        pf->setMarker( 0 );

        // Warning : Vtk orientation is not the same as Feel orientation !

        pf->setPoint( 0, mesh->point( _vtkMesh->GetCell( __i )->GetPointId( 0 ) ) );
        pf->setPoint( 1, mesh->point( _vtkMesh->GetCell( __i )->GetPointId( 1 ) ) );
        pf->setPoint( 2, mesh->point( _vtkMesh->GetCell( __i )->GetPointId( 2 ) ) );
        Debug( 8099 ) << "[FilterFromVtk] point 0 " << pf->point( 0 ).node() << " global id: " << pf->point( 0 ).id() << "\n"
                      << "[FilterFromVtk] point 1 " << pf->point( 1 ).node() << " global id: " << pf->point( 1 ).id() << "\n"
                      << "[FilterFromVtk] point 2 " << pf->point( 2 ).node() << " global id: " << pf->point( 2 ).id() << "\n";
        mesh->addElement( *pf );
        delete pf;
    }


    FEELPP_ASSERT( n_faces == mesh->numFaces() )( n_faces )( mesh->numFaces() ).error( "invalid face container size" );

    Debug( 8099 ) <<"[FilterFromVtk] done with element accumulation !\n";

    mesh->setNumVertices( __n );

    Debug( 8099 ) <<"[FilterFromVtk] Face Update !\n";


    Debug( 8099 ) <<"[FilterFromVtk] Face Update Successful !\n";

#else
    std::cerr << "The library was not compiled with vtk support\n";
#endif /* FEELPP_HAS_VTK */

}

template<typename MeshType>
void
FilterFromVtk3D<MeshType>::visit( mesh_type* mesh, mpl::int_<3> )
{
    detail::ignore_unused_variable_warning( mesh );
#if defined(FEELPP_HAS_VTK)
    //  std::cout <<"Start of mesh conversion !" << std::endl;

    vtkUnstructuredGrid * _vtkMesh = this->getVtkMesh();

    uint __n = _vtkMesh->GetNumberOfPoints(); // Nbre of nodes

    Debug( 8099 ) <<"Number of points : "<< __n << "\n";

    uint __nele = _vtkMesh->GetNumberOfCells(); // Nbre of elements

    Debug( 8099 ) <<"Number of elements : "<< __nele << "\n";

    // add the points to the mesh

    for ( uint __i = 0; __i < __n; ++__i )
    {
        node_type __nd( 3 );
        __nd[0] = _vtkMesh->GetPoint( __i )[0];
        __nd[1] = _vtkMesh->GetPoint( __i )[1];
        __nd[2] = _vtkMesh->GetPoint( __i )[2];
        point_type __pt( __i,__nd, false );

        if ( __nd[0] == -1 || __nd[1] == -1 || __nd[2] == -1 || __nd[0] + __nd[1] + __nd[2] == 0 )

        {
            __pt.setOnBoundary( true );
            __pt.marker() = 0;
        }

        else
        {
            __pt.setOnBoundary( false );
            __pt.marker() = 1;
        }

        mesh->addPoint( __pt );
    }

    Debug( 8099 ) << "[FilterFromVtk3D] mesh np = " << mesh->numPoints() << "\n";
    FEELPP_ASSERT( mesh->numPoints() == __n )( __n )( mesh->numPoints() ).error( "invalid number of points" );
    size_type n_faces = 0;

    // Add Boundary faces

    face_type* pf0 = new face_type;

    pf0->setMarker( 1 );
    pf0->setPoint( 0, mesh->point( 1 ) );
    pf0->setPoint( 1, mesh->point( 3 ) );
    pf0->setPoint( 2, mesh->point( 2 ) );


    pf0->setId( n_faces++ );
    pf0->setOnBoundary( true );
    mesh->addFace( *pf0 );

    delete pf0;

    face_type* pf1 = new face_type;

    pf1->setMarker( 1 );
    pf1->setPoint( 0, mesh->point( 2 ) );
    pf1->setPoint( 1, mesh->point( 3 ) );
    pf1->setPoint( 2, mesh->point( 0 ) );


    pf1->setId( n_faces++ );
    pf1->setOnBoundary( true );
    mesh->addFace( *pf1 );

    delete pf1;

    face_type* pf2 = new face_type;

    pf2->setMarker( 1 );
    pf2->setPoint( 0, mesh->point( 3 ) );
    pf2->setPoint( 1, mesh->point( 1 ) );
    pf2->setPoint( 2, mesh->point( 0 ) );

    pf2->setId( n_faces++ );
    pf2->setOnBoundary( true );
    mesh->addFace( *pf2 );

    delete pf2;

    face_type* pf3 = new face_type;

    pf3->setMarker( 1 );
    pf3->setPoint( 0, mesh->point( 0 ) );
    pf3->setPoint( 1, mesh->point( 1 ) );
    pf3->setPoint( 2, mesh->point( 2 ) );

    pf3->setId( n_faces++ );
    pf3->setOnBoundary( true );
    mesh->addFace( *pf3 );

    delete pf3;

    // add the elements to the mesh

    for ( uint __i = 0; __i < __nele; ++__i )
    {
        // Here we only have triangular elements of order 1

        element_type * pf = new element_type;

        pf->setId( __i );
        pf->setMarker( 0  );

        // Warning : Vtk orientation is not the same as Feel orientation !

        pf->setPoint( 0, mesh->point( _vtkMesh->GetCell( __i )->GetPointId( 0 ) ) );
        pf->setPoint( 1, mesh->point( _vtkMesh->GetCell( __i )->GetPointId( 2 ) ) );
        pf->setPoint( 2, mesh->point( _vtkMesh->GetCell( __i )->GetPointId( 1 ) ) );
        pf->setPoint( 3, mesh->point( _vtkMesh->GetCell( __i )->GetPointId( 3 ) ) );

        mesh->addElement( *pf );
        delete pf;
    }


    FEELPP_ASSERT( n_faces == mesh->numFaces() )( n_faces )( mesh->numFaces() ).error( "invalid face container size" );

    Debug( 8099 ) <<"[FilterFromVtk] done with element accumulation !\n";

    mesh->setNumVertices( __n );

    Debug( 8099 ) <<"[FilterFromVtk] Face Update !\n";

    // do not renumber the mesh entities
    //mesh->updateForUse( MESH_ALL_COMPONENTS & (~MESH_RENUMBER) );

    Debug( 8099 ) <<"[FilterFromVtk] Face Update Successful !\n";

#else
    std::cerr << "The library was not compiled with vtk support\n";
#endif /* FEELPP_HAS_VTK */

}

} // Feel

#endif /* __filter_H */
