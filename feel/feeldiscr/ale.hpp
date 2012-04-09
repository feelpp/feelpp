/*
  This file is part of the Feel library

  Copyright (C) 2007,2008 EPFL
  Copyright (C) 2010 University of Coimbra

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
   \file ale.hpp
   \author Goncalo Pena <goncalo.pena@epfl.ch>
   \date 2008-04-17
*/


#ifndef __ALE
#define __ALE 1



#include <feel/feelcore/feel.hpp>

#include <feel/feelalg/backend.hpp>

#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/interpolate.hpp>
#include <feel/feeldiscr/meshhighorder.hpp>

#include <feel/feelfilters/exporterquick.hpp>

#include <feel/feelmesh/filters.hpp>
#include <feel/feelmesh/meshmover.hpp>

#include <feel/feelpoly/im.hpp>
#include <feel/feelpoly/fekete.hpp>
#include <feel/feelpoly/gausslobatto.hpp>

#include <feel/feelvf/vf.hpp>


/*
 * Class to construct ALE mappings (only works for 2D meshes)
 */

namespace Feel
{
template< class Convex >
class ALE
{
    static const bool is_simplex = Convex::is_simplex;
    static const uint16_type Dim = Convex::nDim;
    static const uint16_type Order = Convex::nOrder;

    /*
     *  Definitions for meshes and functionspaces
     */
    template< int i >
    struct MyConvex
    {
        typedef typename mpl::if_< mpl::bool_< is_simplex >, Simplex<Dim, i>, Hypercube<Dim, i> >::type type;
    };

    typedef typename MyConvex<1>::type convex_type;

    template< int i >
    struct MyPointSet
    {
        typedef PointSetFekete< convex_type, i, double> type;
    };


    typedef Mesh< convex_type > mesh_type;
    typedef boost::shared_ptr<mesh_type> mesh_ptrtype;


    template< int i >
    struct MyBasis
    {
        typedef typename MyPointSet<i>::type pointSet_type;

        typedef bases<Lagrange<i, Vectorial, Continuous, PointSetFekete> > type;
    };



    typedef typename MyBasis<1>::type p1_basis_type;

    typedef FunctionSpace< mesh_type, p1_basis_type, double> p1_functionspace_type;
    typedef boost::shared_ptr<p1_functionspace_type> p1_functionspace_ptrtype;
    typedef typename p1_functionspace_type::element_type p1_element_type;


    typedef typename mesh_type::point_type point_type;
    typedef typename mesh_type::super_faces::marker_face_const_iterator marker_face_const_iterator;
    typedef typename mesh_type::gm_type::points_type points_type;
    typedef typename mesh_type::element_type geo_element_type;

    typedef typename geo_element_type::gm_type gm_type;
    typedef typename geo_element_type::gm_ptrtype gm_ptrtype;
    typedef typename gm_type::template Context<vm::POINT, geo_element_type> gm_context_type;
    typedef boost::shared_ptr<gm_context_type> gm_context_ptrtype;



    /*
     *  Definitions for PN mesh and functionspace
     */
    typedef typename MyConvex<Order>::type new_convex_type;
    typedef Mesh< new_convex_type > new_mesh_type;
    typedef boost::shared_ptr<new_mesh_type> new_mesh_ptrtype;

    typedef typename MyBasis<Order>::type pN_basis_type;

    typedef typename MyPointSet<Order>::type pointset_type;


    /*backend typedefs*/
    typedef Backend<double> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;

    typedef std::pair<double,double> interval_type;

    typedef MeshHighOrder< new_convex_type > ho_mesh_type;
    typedef boost::shared_ptr< ho_mesh_type > ho_mesh_ptrtype;


public:

    typedef FunctionSpace< mesh_type, pN_basis_type, double> pN_functionspace_type;
    typedef boost::shared_ptr<pN_functionspace_type> pN_functionspace_ptrtype;
    typedef typename pN_functionspace_type::element_type pN_element_type;

    typedef typename pN_functionspace_type::fe_type fe_type;
    typedef typename fe_type::template Context<vm::POINT|vm::JACOBIAN, fe_type, gm_type, geo_element_type, gm_context_type::context> fectx_type;
    typedef boost::shared_ptr<fectx_type> fectx_ptrtype;


    /**
     * constructor
     * \param intX pair of doubles describing the inlet and outlet points coordinates (1D)
     */
    ALE( interval_type const& intX, mesh_ptrtype& mesh, po::variables_map const& vm );

    /**
     * copy constructor
     */
    ALE( ALE const& tc );

    /**
     * desctructor
     */
    ~ALE();

    /**
     * \return the high order ALE map
     */
    pN_element_type& getMap();

    /**
     * \return the high order dispacement
     */
    pN_element_type& getDisplacement();

    /**
     * \returns an element containing the position of the points in
     * the reference mesh
     */
    pN_element_type& getIdentity();

    /**
     * \return the reference mesh
     */
    mesh_ptrtype getReferenceMesh();

    new_mesh_ptrtype getMovingMesh();

    pN_functionspace_ptrtype functionSpace();

    template< typename elem_type >
    std::vector<uint16_type> getFacesId( elem_type const&, std::vector<flag_type> const& flagSet );

    points_type getRefFacesPoints( std::vector<uint16_type> const& facesIds,
                                   pointset_type& pointset );

    template< typename elem_type >
    void generateP1BoundaryMap( std::vector<flag_type>& flagSet,
                                std::vector<elem_type> const& referencePolyBoundarySet,
                                std::vector<elem_type> const& polyDisplacementSet,
                                p1_element_type& u );

    template< typename elem_type >
    void generateHighOrderMap( std::vector<flag_type>& flagSet,
                               std::vector<elem_type> const& referencePolyBoundarySet,
                               std::vector<elem_type> const& polyDisplacementSet,
                               bool leaveInternalNodesZero = 0 );

private:

    backend_ptrtype b;

    std::pair<double,double> intervalX;

    mesh_ptrtype reference_mesh;

    p1_functionspace_ptrtype p1_fspace;
    pN_functionspace_ptrtype pN_fspace;

    new_mesh_ptrtype new_mesh;

    sparse_matrix_ptrtype harmonic;

    pN_element_type pN_ale, pN_displacement, pN_identity;

    ho_mesh_ptrtype ho_mesh;

    boost::timer M_timer;

    double M_elast_param;

    void generateP1Map( p1_element_type& p );

    void map2Displacement( pN_element_type const& u, pN_element_type& disp );

    void displacement2Map( pN_element_type const& disp, pN_element_type& u );

    template< typename elem_type >
    void updatePointsInFaces( std::vector<flag_type>& flagSet,
                              std::vector<elem_type> const& referencePolyBoundarySet,
                              std::vector<elem_type> const& polyDisplacementSet,
                              pN_element_type& p, mpl::bool_<true> );

    template< typename elem_type >
    void updatePointsInFaces( std::vector<flag_type>& flagSet,
                              std::vector<elem_type> const& referencePolyBoundarySet,
                              std::vector<elem_type> const& polyDisplacementSet,
                              pN_element_type& p, mpl::bool_<false> )
    {
    }


};

template < class Convex >
template< typename elem_type >
std::vector<uint16_type>
ALE<Convex>::getFacesId( elem_type const& elt, std::vector<flag_type> const& flagSet )
{
    std::vector<uint16_type> facesId;

    for ( uint8_type i = 0; i < geo_element_type::numEdges; ++i )
    {
        auto edge = elt.edge( i );

        auto result = find( flagSet.begin(),
                            flagSet.end(),
                            edge.marker().value() );

        if ( result != flagSet.end() )
        {
#if !defined ( NDEBUG )
            Debug( 1234 ) << "We are in element " << elt.id() << "\n";
            Debug( 1234 ) << "Marker of edge : " << edge.marker().value() << "\n";
#endif
            facesId.push_back( i );
        }
    }

    return facesId;
}

template < class Convex >
typename ALE<Convex>::points_type
ALE<Convex>::getRefFacesPoints( std::vector<uint16_type> const& facesId,
                                pointset_type& pointset )
{
    // generate points in reference element and apply Gordon Hall transformation
    points_type pts ( Dim, facesId.size()*pointset_type::nbPtsPerEdge + pointset_type::nbPtsPerFace );
    auto it_faces = facesId.begin();
    auto en_faces = facesId.end();


    uint16_type counter = 0;

    for ( ; it_faces != en_faces; ++it_faces )
    {
        ublas::subrange( pts, 0, Dim, counter, counter+pointset_type::nbPtsPerEdge ) =
            pointset.pointsBySubEntity( Dim-1, *it_faces );

        counter += pointset_type::nbPtsPerEdge;
    }

    ublas::subrange( pts, 0, Dim, counter, pts.size2() ) = pointset.pointsBySubEntity( Dim, 0 );
    return pts;
}
template < class Convex >
template< typename elem_type >
void
ALE<Convex>::generateHighOrderMap( std::vector<flag_type>& flagSet,
                                   std::vector<elem_type> const& referencePolyBoundarySet,
                                   std::vector<elem_type> const& polyDisplacementSet,
                                   bool leaveInternalNodesZero )
{
    using namespace Feel::vf;

    M_timer.restart();
    // first, we generate the p1 ale map
    p1_element_type p1_displacement( p1_fspace, "p1_displacement" );
    generateP1BoundaryMap( flagSet, referencePolyBoundarySet, polyDisplacementSet, p1_displacement );
    //Log() << "[ALE] Time to generate P1 boundary map: " << M_timer.elapsed() << "\n";

    if ( !leaveInternalNodesZero )
    {
        M_timer.restart();
        generateP1Map( p1_displacement );
        //Log() << "[ALE] Time to generate P1 map: " << M_timer.elapsed() << "\n";
    }

    M_timer.restart();
    interpolate( pN_fspace, p1_displacement, pN_displacement, INTERPOLATE_SAME_MESH );
    //Log() << "[ALE] Time to generate PN projection of P1 map: " << M_timer.elapsed() << "\n";

#if !defined ( NDEBUG )
    MeshHighOrder< convex_type > auxiliar_mesh ( reference_mesh );
    auxiliar_mesh.generateMesh( flagSet, referencePolyBoundarySet );
    mesh_ptrtype aux_mesh = auxiliar_mesh.getMesh();

    MeshMover<mesh_type> p1_mesh_mover;
    p1_mesh_mover.apply( aux_mesh, p1_displacement );

    ExporterQuick<mesh_type> exp( "ALE", "ensight" );
    exp.save( aux_mesh );

    MeshHighOrder< convex_type > auxiliar_mesh2 ( reference_mesh );
    auxiliar_mesh2.generateMesh( flagSet, referencePolyBoundarySet );
    mesh_ptrtype aux_mesh2 = auxiliar_mesh2.getMesh();

    p1_element_type p1_ex_displacement( p1_fspace, "p1_ex_displacement" );
    p1_element_type p1_ex_diff( p1_fspace, "p1_ex_diff" );
    AUTO( f,  ( 0.02096834639998*Px()*Px()*Px()    -0.09351398506873*Px()*Px()    -0.09961900352798*Px()+    1.3000 )  ) ;
    AUTO( f2, ( -0.0210*Px()*Px()*Px() + 0.0935*Px()*Px() + 0.0996*Px() -1.4000 )  );
    p1_ex_displacement = vf::project( p1_fspace, markedfaces( p1_fspace->mesh(),3 ), vec( Px(),f ) );
    p1_ex_diff = vf::project( p1_fspace, markedfaces( p1_fspace->mesh(),3 ), idv( pN_displacement )-vec( constant( 0. ),f ) );
    //p1_ex_displacement = vf::project( p1_fspace, markedfaces(mesh,1), vec(Px(),f2) );
    //exp.save( 0, p1_displacement, pN_displacement );

    //MeshMover<mesh_type> p1_mesh_mover2;
    //p1_mesh_mover2.apply(aux_mesh2, p1_ex_displacement);

    ExporterQuick<mesh_type> exp_ex( "ALE_ex", "ensight" );
    exp_ex.save( 0, p1_ex_displacement, p1_ex_diff );

#endif

    displacement2Map( pN_displacement, pN_ale );
    updatePointsInFaces( flagSet, referencePolyBoundarySet, polyDisplacementSet, pN_ale, mpl::bool_< ( Order > 1 ) >() );
    map2Displacement( pN_ale, pN_displacement );

}


template < class Convex >
template< typename elem_type >
void
ALE<Convex>::generateP1BoundaryMap( std::vector<flag_type>& flagSet,
                                    std::vector<elem_type> const& referencePolyBoundarySet,
                                    std::vector<elem_type> const& polyDisplacementSet,
                                    p1_element_type& u )
{
    using namespace Feel::vf;

    double tolerance = 1e-10;

    double Y_coordinate_first = 0;
    double Y_coordinate_second = 0;


    // define boundary conditions
    // first, we start with the part of the boundary that is fixed
    typedef typename elem_type::functionspace_type::node_type node_type;

    node_type pt( 1 );
    pt[0] = intervalX.first;

    elem_type p1 = polyDisplacementSet[0];
    elem_type p2 = polyDisplacementSet[1];

    elem_type p1_aux = referencePolyBoundarySet[0];
    elem_type p2_aux = referencePolyBoundarySet[1];

    double pt1_aux = p1_aux( pt )( 0,0,0 );
    double pt2_aux = p2_aux( pt )( 0,0,0 );

    double pt1 = p1( pt )( 0,0,0 );
    double pt2 = p2( pt )( 0,0,0 );

    Y_coordinate_first = pt1_aux;
    Y_coordinate_second = pt2_aux;

    double m = ( pt2-pt1 )/( Y_coordinate_second - Y_coordinate_first );

    AUTO( aux, ( cst_ref( pt1 ) + cst_ref( m )*( Py() - Y_coordinate_first ) )*oneY() );

    u = vf::project( p1_fspace, boundaryelements( *reference_mesh ), aux*chi( Px() < intervalX.first + tolerance ) );

    pt[0] = intervalX.second;

    pt1 = p1( pt )( 0,0,0 );
    pt2 = p2( pt )( 0,0,0 );
    pt1_aux = p1_aux( pt )( 0,0,0 );
    pt2_aux = p2_aux( pt )( 0,0,0 );

    Y_coordinate_first = pt1_aux;
    Y_coordinate_second = pt2_aux;

    m = ( pt2-pt1 )/( Y_coordinate_second - Y_coordinate_first );

    u += vf::project( p1_fspace, boundaryelements( *reference_mesh ), aux*chi( Px() > intervalX.second - tolerance ) );

    std::map< flag_type, std::vector<uint16_type> > pointIdOnBoundary;

    double pointNode0 = 0;
    double pointNode1 = 0;
    double polyNode = 0;

    AUTO( f1, Px() - cst_ref( pointNode0 ) );
    AUTO( f2, Py() - cst_ref( pointNode1 ) );
    AUTO( radius, f1*f1 + f2*f2 );
    AUTO( dirac_function, cst_ref( polyNode )*oneY()*chi( radius < tolerance*tolerance ) );

    for ( uint16_type i=0; i < flagSet.size(); ++i )
    {
        marker_face_const_iterator face_it = reference_mesh->beginFaceWithMarker( flagSet[i] );
        marker_face_const_iterator face_en = reference_mesh->endFaceWithMarker( flagSet[i] );

        for ( ; face_it != face_en; ++face_it )
        {
            for ( uint16_type j=0; j < 2; ++j )
            {
                point_type point = reference_mesh->face( face_it->id() ).point( j );

                std::vector<uint16_type>::const_iterator result;

                result = find( pointIdOnBoundary[ flagSet[i] ].begin(),
                               pointIdOnBoundary[ flagSet[i] ].end(), point.id() );

                if ( ( result == pointIdOnBoundary[ flagSet[i] ].end() ) || pointIdOnBoundary[ flagSet[i] ].empty() )
                {
                    ublas::vector<double> pt_coord ( 2 );
                    pt_coord[0] = ( point.node() )[0];

                    if ( ( pt_coord[0] > intervalX.first + tolerance ) && ( pt_coord[0] < intervalX.second - tolerance ) )
                    {
                        p1 = polyDisplacementSet[i];

                        pt[0] = ( point.node() )[0];
                        polyNode = p1( pt )( 0,0,0 );

                        pointNode0 = pt[0];
                        pointNode1 = point.node()[1];

#if !defined ( NDEBUG )
                        Debug( 1234 ) << "From (" << pt[0] << "," << ( point.node() )[1]
                                      << ") to (" << pt[0] << "," << polyNode << ")\n";
#endif

                        u += vf::project( p1_fspace,
                                          idedelements( *reference_mesh, face_it->ad_first()  ),
                                          dirac_function );

                        pointIdOnBoundary[ flagSet[i] ].push_back( point.id() );
                    }
                }
            }
        }
    }
}


template < class Convex >
template< typename elem_type >
void
ALE<Convex>::updatePointsInFaces( std::vector<flag_type>& flagSet,
                                  std::vector<elem_type> const& referencePolyBoundarySet,
                                  std::vector<elem_type> const& polyDisplacementSet,
                                  pN_element_type& p, mpl::bool_<true> )
{
    M_timer.restart();

    using namespace Feel::vf;

    pointset_type pointset;

    double tolerance = 1e-10;

    typename pN_functionspace_type::node_type node( Dim );
    AUTO( f1, Px() - cst_ref( node[0] ) );
    AUTO( f2, Py() - cst_ref( node[1] ) );
    AUTO( radius, f1*f1 + f2*f2 );


    gm_ptrtype gm( new gm_type );
    typename gm_type::precompute_ptrtype __geopc( new typename gm_type::precompute_type( gm, pointset.points() ) );

    gm_context_ptrtype __c( new gm_context_type( gm, reference_mesh->element( 0 ), __geopc ) );

    boost::timer time;

    std::vector<elem_type> polyBoundarySet;

    for ( uint16_type i = 0; i < referencePolyBoundarySet.size(); ++i )
    {
        elem_type temp ( referencePolyBoundarySet[i].functionSpace(), "temp" );
        temp = referencePolyBoundarySet[i];
        temp += polyDisplacementSet[i];

        polyBoundarySet.push_back( temp );
    }

    auto it_elt = reference_mesh->beginElementOnBoundary();
    auto en_elt = reference_mesh->endElementOnBoundary();
    typedef boost::multi_array<double,3> array_type;
    auto facesId = this->getFacesId( *it_elt, flagSet );
    typedef typename pN_element_type::pc_type pN_pc_type;
    typedef typename pN_element_type::pc_ptrtype pN_pc_ptrtype;
    pN_pc_ptrtype pc( new pN_pc_type( pN_ale.functionSpace()->fe(), pointset.points() ) );
    fectx_ptrtype fectx( new fectx_type( pN_ale.functionSpace()->fe(),
                                         __c,
                                         pc ) );
    array_type interpfunc( pN_ale.idExtents( *fectx ) );
    array_type interpfunc_vertices( pN_ale.idExtents( *fectx ) );

    for ( ; it_elt != en_elt; ++it_elt )
    {
        // find if element has a face in the boundary
        std::vector<uint16_type> facesId = this->getFacesId( *it_elt, flagSet );

        if ( facesId.size() > 0  )
        {
            // generate points in reference element
            auto pts = getRefFacesPoints( facesId, pointset );

            // generate points in real element
            __geopc->update( pts );
            __c->update( reference_mesh->element( it_elt->id() ), __geopc );
            pc->update( __c->xRefs() );
            fectx->update( __c, pc );
            std::fill( interpfunc.data(),
                       interpfunc.data()+interpfunc.num_elements(), 0 );
            pN_ale.id( *fectx, interpfunc );
            points_type pts_reference = __c->xReal();

            __geopc->update( pointset.pointsByEntity( 0 ) );
            __c->update( reference_mesh->element( it_elt->id() ), __geopc );
            pc->update( __c->xRefs() );
            fectx->update( __c, pc );
            std::fill( interpfunc_vertices.data(),
                       interpfunc_vertices.data()+interpfunc_vertices.num_elements(), 0 );
            pN_ale.id( *fectx, interpfunc_vertices );

            geo_element_type copy_elt = *it_elt;
            ublas::vector<double> v( 2 );
            std::vector< ublas::vector<double> > vertices( geo_element_type::numVertices );

            for ( uint16_type i=0; i < geo_element_type::numVertices; ++i )
            {
                v( 0 ) = interpfunc_vertices[0][0][i] - ( it_elt->point( i ).node() )[0];
                v( 1 ) = interpfunc_vertices[1][0][i] - ( it_elt->point( i ).node() )[1];

                copy_elt.applyDisplacement( i, v );
                vertices[i] = v;
            }

            // apply Gordon Hall transformation
            ho_mesh->GordonHall( *it_elt, pts, flagSet, polyBoundarySet );

            for ( uint16_type i=0; i < geo_element_type::numVertices; ++i )
            {
                v( 0 ) = -( interpfunc_vertices[0][0][i] - ( it_elt->point( i ).node() )[0] );
                v( 1 ) = -( interpfunc_vertices[1][0][i] - ( it_elt->point( i ).node() )[1] );
                v = -vertices[i];
                copy_elt.applyDisplacement( i, v );
            }

            for ( uint16_type j = 0; j < pts.size2(); ++j )
            {
                node[0] = pts_reference( 0,j );
                node[1] = pts_reference( 1,j );

#if !defined ( NDEBUG )
                Debug( 1234 ) << "Point (" << node[0] << "," << node[1]
                              << ") moves to (" << pts( 0,j ) << "," << pts( 1,j )
                              << ")" << "\n";
#endif

                p += vf::project( pN_fspace, idedelements( *reference_mesh, it_elt->id() ),
                                  ( ( pts( 0,j ) - interpfunc[0][0][j] )*oneX() +
                                    ( pts( 1,j ) - interpfunc[1][0][j] )*oneY() )*
                                  chi( ( radius < tolerance*tolerance ) ) );
            }
        }
    }

    Log() << "[ALE] Time to update PN map (faces): " << M_timer.elapsed() << "\n";
}

} // Feel

#if !defined( FEELPP_INSTANTIATION_MODE )
# include <feel/feeldiscr/ale.cpp>
#endif // FEELPP_INSTANTIATION_MODE
#endif // __ALE
