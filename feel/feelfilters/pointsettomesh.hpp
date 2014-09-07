/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-11-25

  Copyright (C) 2005,2006 EPFL

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
   \file pointsettomesh.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-11-25
 */
#ifndef __PointSetToMesh_H
#define __PointSetToMesh_H 1

#include <boost/shared_ptr.hpp>
#include <boost/optional.hpp>

#include <feel/feelcore/feel.hpp>

#include <stdlib.h>

#if defined(FEELPP_HAS_VTK)

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-W#warnings"
#endif

#include <vtkVersion.h>
#include <vtkPointSet.h>
#include <vtkDelaunay2D.h>
#include <vtkDelaunay3D.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkExtractEdges.h>
#include <vtkPolyDataMapper.h>
#include <vtkDataSetMapper.h>

#if defined(__clang__)
#pragma clang diagnostic pop
#endif

#endif /* FEELPP_HAS_VTK */

#include <feel/feelmesh/pointset.hpp>
#include <feel/feelfilters/filterfromvtk.hpp>

namespace Feel
{

/**
 * \class PointSetToMesh
 * \brief transform a point set to a mesh data structure using a Delaunay
 *
 * the delaunay algorithm comes from Vtk
 *
 * @author Christophe Prud'homme
 * @see
 */
template<typename Convex, typename T>
class PointSetToMesh
    :
public VisitorBase,
public Visitor<PointSet<Convex, T> >
{
    typedef VisitorBase super1;
    typedef Visitor<PointSet<Convex, T> > super2;

public:
    /** @name Constants
     */
    //@{

    static const uint16_type nDim = Convex::nDim;

    //@}

    /** @name Typedefs
     */
    //@{

    typedef T value_type;
    typedef Convex convex_type;
    typedef typename convex_type::template shape<convex_type::nDim, 1, convex_type::nDim>::type mesh_convex_type;

    typedef PointSet<convex_type, value_type> pointset_type;

    typedef Mesh<mesh_convex_type> mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;

    typedef typename matrix_node<value_type>::type points_type;

    typedef typename mesh_type::point_type point_type;
    typedef typename point_type::node_type node_type;
    typedef typename mesh_type::edge_type edge_type;
    typedef typename mesh_type::face_type face_type;
    typedef typename mesh_type::element_type element_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    PointSetToMesh()
        :
        super1(),
        super2(),
        M_mesh( new mesh_type( Environment::worldCommSeq() ) ),
        M_vertices()
    {}
    PointSetToMesh( PointSetToMesh const & p )
        :
        super1( p ),
        super2( p ),
        M_mesh( p.M_mesh ),
        M_vertices( p.M_vertices )
    {}

    ~PointSetToMesh()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    mesh_ptrtype mesh()
    {
        return M_mesh;
    }

    //@}

    /** @name  Mutators
     */
    //@{

    void addBoundaryPoints( points_type const& v )
    {
        M_vertices = v;
    }

    //@}

    /** @name  Methods
     */
    //@{

    /**
     * visit the point set \p pset and transform it into a mesh data
     * structure using Delaunay triangulation
     */
    void visit( pointset_type* pset )
    {
        visit( pset, mpl::int_<nDim>() );
    }


    //@}



protected:

private:

    void visit( pointset_type* pset, mpl::int_<1> );
    void visit( pointset_type* pset, mpl::int_<2> );
    void visit( pointset_type* pset, mpl::int_<3> );
private:

    mesh_ptrtype M_mesh;
    boost::optional<points_type> M_vertices;
};

template<typename Convex, typename T>
void
PointSetToMesh<Convex, T>::visit( pointset_type* pset, mpl::int_<1> )
{
    DVLOG(2) << "[PointSetToMesh::visit(<1>)] pointset to mesh\n";
    M_mesh = mesh_ptrtype( new mesh_type( Environment::worldComm().subWorldCommSeq() ) );

    size_type __npts = pset->nPoints();

    for ( uint __i = 0; __i < __npts; ++__i )
    {
        node_type __nd( 1 );
        __nd = pset->point( __i );

        point_type __pt( __i,__nd, false );
        __pt.marker() = 0;

        if ( __nd[0] == -1 || __nd[0] == 1 )
        {
            __pt.setOnBoundary( true );
        }

        else
        {
            __pt.setOnBoundary( false );

        }

        M_mesh->addPoint( __pt );
    }

    size_type n_faces = 0;


    // Add Boundary faces

    face_type* pf0 = new face_type;

    pf0->marker() = 0;
    pf0->setPoint( 0, M_mesh->point( 0 ) );

    pf0->setId( n_faces++ );
    pf0->setOnBoundary( true );
    M_mesh->addFace( *pf0 );

    delete pf0;

    face_type* pf1 = new face_type;

    pf1->marker() = 0;
    pf1->setPoint( 0, M_mesh->point( 1 ) );

    pf1->setId( n_faces++ );
    pf1->setOnBoundary( true );
    M_mesh->addFace( *pf1 );

    delete pf1;

    size_type __nele = __npts-1;

    for ( uint __i = 0; __i < __nele; ++__i )
    {
        element_type  pf;

        pf.setId( __i );
        pf.marker() = 0;

        if ( __nele == 1 )
        {
            pf.setPoint( 0, M_mesh->point( __i ) );
            pf.setPoint( 1, M_mesh->point( __i+1 ) );
        }

        else if ( __i == 0 )
        {
            pf.setPoint( 0, M_mesh->point( __i ) );
            pf.setPoint( 1, M_mesh->point( __i+2 ) );
        }

        else if ( __i == __nele-1 )
        {
            pf.setPoint( 0, M_mesh->point( __i+1 ) );
            pf.setPoint( 1, M_mesh->point( 1 ) );
        }

        else
        {
            pf.setPoint( 0, M_mesh->point( __i+1 ) );
            pf.setPoint( 1, M_mesh->point( __i+2 ) );
        }

        element_type const& e  = M_mesh->addElement( pf );
        DVLOG(2) << "o element " << e.id() << "\n"
                      << "  p1 = " << e.point( 0 ).node() << "\n"
                      << "  p2 = " << e.point( 1 ).node() << "\n";
    }


    FEELPP_ASSERT( n_faces == M_mesh->numFaces() )( n_faces )( M_mesh->numFaces() ).error( "invalid face container size" );

    DVLOG(2) <<"[PointSetToMesh<1>] done with element accumulation !\n";

    M_mesh->setNumVertices( __npts );

    DVLOG(2) <<"[PointSetToMesh<1>] Face Update !\n";

    //// do not renumber the M_mesh entities
    //M_mesh->updateForUse( MESH_ALL_COMPONENTS & (~MESH_RENUMBER) );

    DVLOG(2) <<"[PointSetToMesh<1>] Face Update Successful !\n";


    DVLOG(2) << "[PointSetToMesh::visit(<1>)] done\n";
}
template<typename Convex, typename T>
void
PointSetToMesh<Convex, T>::visit( pointset_type* pset, mpl::int_<2> )
{
#if defined(FEELPP_HAS_VTK)
    // reinitialize mesh
    M_mesh = mesh_ptrtype( new mesh_type( Environment::worldComm().subWorldCommSeq() ) );

    vtkPoints *newPoints = vtkPoints::New();

    if ( M_vertices )
    {
        DVLOG(2) << "adding vertices\n" << M_vertices.get() << "\n";

        for ( size_type i = 0; i < M_vertices->size2(); ++i )
        {
            newPoints->InsertNextPoint( M_vertices.get()( 0, i ), M_vertices.get()( 1, i ), 0 );
        }
    }

    std::vector<double> perturbation( pset->nPoints() );

    for ( size_type i=0 ; i< pset->nPoints() ; ++i )
    {
        perturbation[i] = std::rand()*1e-6/RAND_MAX;
        uint16_type index = newPoints->InsertNextPoint( pset->points()( 0,i )+perturbation[i], pset->points()( 1,i ), 0 );
        DVLOG(2) << "Inserting point with id " << index << "\n";
        DVLOG(2) << "pset.point( " << i << " )= " << pset->point( i ) << "\n";
    }

    vtkPolyData *polyData = vtkPolyData::New();
    polyData->SetPoints( newPoints );

    vtkDelaunay2D *delaunay2D = vtkDelaunay2D::New();

#if VTK_MAJOR_VERSION <= 5
    delaunay2D->SetInput( polyData );
#else
    delaunay2D->SetInputData( polyData );
#endif

    /**
     * The Offset parameter helps to get a convex hull of the points.
     * If the result mesh is not convex, increase the Offset !
     **/

    //delaunay2D->SetOffset( 8 ) ;
    delaunay2D->Update();

    vtkPolyData* outMesh = delaunay2D->GetOutput( );

    /** General informations about vtkPolyData **/

    outMesh->BuildLinks();

    int Nelem = outMesh->GetNumberOfPolys();
    int Npts = outMesh->GetNumberOfPoints();

    DVLOG(2) << "Number of cells  = " << Nelem   << "\n";
    DVLOG(2) << "Number of points  = " << Npts   << "\n";

    for ( int i=0; i< Npts; ++i )
    {
        // revert back the perturbation
        outMesh->GetPoints()->SetPoint( i,
                                        outMesh->GetPoints()->GetPoint( i )[0]- perturbation[i],
                                        outMesh->GetPoints()->GetPoint( i )[1],
                                        outMesh->GetPoints()->GetPoint( i )[2] );
    }

    for ( int i=0; i< Nelem; ++i )
    {
        //  std::cout << "\nLa cellule num�ro : " << i << " compte " << outMesh->GetCell(i)->GetNumberOfPoints() << " points." << "\n";
        DVLOG(2) << "Element Id = " << i << "\n";
        DVLOG(2) << "Point 0 (" <<  ( int )outMesh->GetCell( i )->GetPointId( 0 ) <<") =" ;
        DVLOG(2) << "(" << outMesh->GetCell( i )->GetPoints()->GetPoint( 0 )[0] << " , "
                      << outMesh->GetCell( i )->GetPoints()->GetPoint( 0 )[1]<< ")" << "\n";
        DVLOG(2) << "Point 1 (" <<  ( int )outMesh->GetCell( i )->GetPointId( 1 ) <<") =" ;
        DVLOG(2) << "(" << outMesh->GetCell( i )->GetPoints()->GetPoint( 1 )[0] << " , "
                      << outMesh->GetCell( i )->GetPoints()->GetPoint( 1 )[1]<< ")" << "\n";
        DVLOG(2) << "Point 2 (" <<  ( int )outMesh->GetCell( i )->GetPointId( 2 ) <<") =" ;
        DVLOG(2) << "(" << outMesh->GetCell( i )->GetPoints()->GetPoint( 2 )[0] << " , "
                      << outMesh->GetCell( i )->GetPoints()->GetPoint( 2 )[1]<< ")" << "\n";

        DVLOG(2) << outMesh->GetCell( i )->GetNumberOfEdges() << "\n";
        DVLOG(2) << ( int )outMesh->GetCell( i )->GetEdge( 0 )->GetPointId( 0 ) << "\n";
        DVLOG(2) << ( int )outMesh->GetCell( i )->GetEdge( 0 )->GetPointId( 1 ) << "\n";
    }



    DVLOG(2) << "[PointSetToMesh::visit(<2>)] delaunay done, now vtk to Mesh<>\n";
    FilterFromVtk<mesh_type> meshfromvtk( outMesh );
    meshfromvtk.visit( M_mesh.get() );
    DVLOG(2) << "[PointSetToMesh::visit(<2>)] done\n";

#else
    std::cerr << "The library was not compiled with vtk support\n";
#endif /* FEELPP_HAS_VTK */

}

template<typename Convex, typename T>
void
PointSetToMesh<Convex, T>::visit( pointset_type* pset, mpl::int_<3> )
{
#if defined(FEELPP_HAS_VTK)
    // reinitialize mesh
    M_mesh = mesh_ptrtype( new mesh_type( Environment::worldComm().subWorldCommSeq() ) );

    vtkPoints *newPoints = vtkPoints::New();

    //if ( M_add_bdy_pts )
    if ( 0 )
    {
        newPoints->InsertNextPoint( -1, -1, -1 );
        newPoints->InsertNextPoint(  1, -1, -1 );
        newPoints->InsertNextPoint( -1,  1, -1 );
        newPoints->InsertNextPoint( -1, -1,  1 );
    }

    for ( size_type i=0 ; i< pset->nPoints() ; ++i )
    {
        newPoints->InsertNextPoint( pset->points()( 0,i ), pset->points()( 1,i ), pset->points()( 2,i ) );
    }

    // create more points

    vtkPolyData *polyData = vtkPolyData::New();
    polyData->SetPoints( newPoints );

    vtkDelaunay3D *delaunay3D = vtkDelaunay3D::New();
#if VTK_MAJOR_VERSION <= 5
    delaunay3D->SetInput( polyData );
#else
    delaunay3D->SetInputData( polyData );
#endif
    delaunay3D->SetOffset( 5 );

    DVLOG(2) <<"[PointSetToMesh::visit(<3>)] Offset = " << delaunay3D->GetOffset() << "\n";
    delaunay3D->Update();

    vtkUnstructuredGrid* outMesh = delaunay3D->GetOutput( );


    DVLOG(2) << "[PointSetToMesh::visit(<3>)] delaunay done, now vtk to Mesh<>\n";
    FilterFromVtk3D<mesh_type> meshfromvtk( outMesh );
    meshfromvtk.visit( M_mesh.get() );
    DVLOG(2) << "[PointSetToMesh::visit(<3>)] done\n";


#else
    std::cerr << "The library was not compiled with vtk support\n";
#endif /* FEELPP_HAS_VTK */
}
} // Feel
#endif /* __PointSetToMesh_H */
