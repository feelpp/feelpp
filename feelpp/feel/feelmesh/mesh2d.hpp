// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 05 Feb 2017
//! @copyright 2017 Feel++ Consortium
//!
#ifndef FEELPP_MESH2D_HPP
#define FEELPP_MESH2D_HPP 1

#include <cstdlib>
#include <fstream>
#include <iomanip>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <boost/foreach.hpp>
#include <boost/mpl/print.hpp>
#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/shared_ptr.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/visitor.hpp>

#include <feel/feelmesh/meshbase.hpp>

#include <feel/feelmesh/geoelement.hpp>

#include <feel/feelmesh/elements.hpp>
#include <feel/feelmesh/faces.hpp>
#include <feel/feelmesh/functors.hpp>
#include <feel/feelmesh/points.hpp>

namespace Feel
{
//!
//! @brief 2D mesh class
//! @ingroup Mesh
//!
//! \code
//! // create a 2D mesh made of simplex of order 1
//! Mesh2D<Simplex<2,1> > mesh;
//!
//! // create a 2D mesh made of simplex of order 2
//! Mesh2D<Simplex<2,2> > mesh;
//! \endcode
//!
//! @author Christophe Prud'homme
//!  @see
//!
template <typename Shape, typename T = double, typename IndexT = uint32_type>
class Mesh2D
    : public VisitableBase<>,
      public MeshBase<IndexT>,
      public Elements<Shape,T>,
      public Points<Shape::nRealDim,T>,
      public Faces<typename Shape::template shape<1, Shape::nOrder, Shape::nRealDim>,
                   typename Elements<Shape,T>::element_type>
{
    // check at compilation time that the shape has indeed dimension 2
    BOOST_STATIC_ASSERT( Shape::nDim == 2 );

  public:
    /** @name Typedefs
     */
    //@{

    static const uint16_type nDim = Shape::nRealDim;
    static const uint16_type nRealDim = Shape::nRealDim;

    typedef typename VisitableBase<>::return_type return_type;

    typedef VisitableBase<> super_visitable;
    typedef MeshBase<IndexT> super;

    using index_type = typename super::index_type;
    using size_type = typename super::size_type;

    typedef Elements<Shape,T> super_elements;
    typedef typename super_elements::elements_type elements_type;
    typedef typename super_elements::element_type element_type;
    typedef typename super_elements::element_iterator element_iterator;
    typedef typename super_elements::element_const_iterator element_const_iterator;
    typedef typename super_elements::update_element_neighbor_type update_element_neighbor_type;

    typedef Points<Shape::nRealDim,T> super_points;
    typedef typename super_points::points_type points_type;
    typedef typename super_points::point_type point_type;

    typedef Faces<typename Shape::template shape<1, Shape::nOrder, Shape::nRealDim>,
                  typename super_elements::element_type>
        super_faces;
    typedef typename super_faces::face_iterator face_iterator;
    typedef typename super_faces::faces_type faces_type;
    typedef typename super_faces::face_type face_type;
    typedef face_type edge_type;

    typedef super_faces super_edges;
    typedef typename super_edges::face_iterator edge_iterator;
    typedef typename super_edges::face_const_iterator edge_const_iterator;
    typedef typename super_edges::faces_reference_wrapper_type edges_reference_wrapper_type;
    typedef typename super_edges::faces_reference_wrapper_ptrtype edges_reference_wrapper_ptrtype;
    typedef typename edges_reference_wrapper_type::iterator edge_reference_wrapper_iterator;
    typedef typename edges_reference_wrapper_type::const_iterator edge_reference_wrapper_const_iterator;

    typedef Mesh2D<Shape,T> self_type;
    typedef std::shared_ptr<self_type> self_ptrtype;

    typedef typename element_type::edge_permutation_type edge_permutation_type;
    typedef typename element_type::face_permutation_type face_permutation_type;

    typedef typename super::face_processor_type face_processor_type;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    /**
     * default constructor
     */
    explicit Mesh2D( std::string const& name = "", worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() )
        : super_visitable(),
          super( name, 2, nRealDim, worldComm ),
          super_elements( worldComm ),
          super_points( worldComm ),
          super_faces( worldComm )
    {
        DVLOG( 2 ) << "[Mesh2D] constructor...\n";
    }

    Mesh2D( Mesh2D const& m ) = default;
    Mesh2D( Mesh2D&& m ) = default;

    /**
     * destructor
     */
    ~Mesh2D() override
    {
        VLOG( 1 ) << "Mesh2D destructor";
        this->clear();
    }

    //@}

    /** @name Operator overloads
     */
    //@{

    Mesh2D& operator=( Mesh2D const& m ) = default;
    Mesh2D& operator=( Mesh2D&& m ) = default;

    //@}

    /** @name Accessors
     */
    //@{

    /**
     * \return \p true if all containers are empty, \p false otherwise
     */
    bool isEmpty() const override
    {
        return ( super_elements::isEmpty() &&
                 super_points::isEmpty() &&
                 super_faces::isEmpty() );
    }

    /**
     * \return the number of elements
     */
    size_type numElements() const override
    {
        return this->elements().size();
    }

    /**
 * \return the number of faces in an element
 * \sa numLocalFaces()
 */
    size_type numLocalEdges() const
    {
        return super_elements::element_type::numLocalEdges;
    }

    /**
     * \return the number of faces in an element
     * \sa numLocalEdges()
     */
    size_type numLocalFaces() const
    {
        return super_elements::element_type::numLocalEdges;
    }

    //!
    //! the number of topological faces per element
    //! @return the number of topological faces per element
    //!
    uint16_type numLocalTopologicalFaces() const
    {
        return super_elements::element_type::numTopologicalFaces;
    }

    /**
     * \return the number of vertices in an element
     */
    size_type numLocalVertices() const
    {
        return super_elements::element_type::numLocalVertices;
    }

    /**
     * \return the number of faces
     */
    size_type numFaces() const override
    {
        return this->faces().size();
    }

    /**
 * \return the number of edges
 */
    size_type numEdges() const
    {
        return this->faces().size();
    }

    /**
 * \return the number of points
 */
    size_type numPoints() const override
    {
        return this->points().size();
    }

    //@}

    /** @name  Mutators
 */
    //@{

    //@}

    /** @name  Methods
 */
    //@{

    face_iterator beginEdge() { return this->beginFace(); }
    face_iterator endEdge() { return this->endFace(); }
    faces_type edges() { return this->faces(); }

    void setWorldComm( worldcomm_ptr_t const& _worldComm ) override
    {
        MeshBase<IndexT>::setWorldComm( _worldComm );
        this->setWorldCommElements( _worldComm );
        this->setWorldCommFaces( _worldComm );
        this->setWorldCommPoints( _worldComm );
    }

    /**
 * clear out all data from the mesh, \p isEmpty() should return
 * \p true after a \p clear()
 */
    void clear() override
    {
        VLOG( 1 ) << "Deleting Mesh2D...\n";

        super::clear();
        super_elements::clear();
        super_points::clear();
        super_faces::clear();
        CHECK( isEmpty() ) << "all mesh containers should be empty after a clear.";
    }

    FEELPP_DEFINE_VISITABLE();
    //@}

  protected:
    /**
     * dummy  implementation
     * \see Mesh
     */
    void renumber() override
    {
        FEELPP_ASSERT( 0 )
            .error( "invalid call" );
    }

    /**
     * update permutation of entities of co-dimension 1
     */
    void updateEntitiesCoDimensionOnePermutation()
        {
            std::vector<index_type> pointIdInElt0(element_type::numVertices);
            std::vector<index_type> pointIdInElt1(element_type::numVertices);

            for ( typename super_elements::element_iterator elt_it = this->beginElement();
                  elt_it != this->endElement(); ++elt_it )
            {
                auto & elt = elt_it->second;
                for ( uint16_type f = 0; f < element_type::numVertices; ++f )
                    pointIdInElt0[f] = elt.point( f ).id();

                for ( uint16_type j = 0; j < element_type::numEdges; j++ )
                {
                    if ( !elt.hasFace( j ) )
                        continue;
                    auto const& face = elt.face( j );

                    if ( !(face.isConnectedTo0() && face.isConnectedTo1()) )
                        continue;

                    bool eltUseConnection0 = face.element0().id() == elt.id();
                    // only change in connection 1
                    if ( eltUseConnection0 )
                        continue;

                    uint16_type idInElement0 = eltUseConnection0? face.idInElement0() : face.idInElement1();
                    uint16_type idInElement1 = eltUseConnection0? face.idInElement1() : face.idInElement0();
                    auto const& elt0 = elt;
                    auto const& elt1 = eltUseConnection0? face.element1() : face.element0();

                    for ( uint16_type f = 0; f < element_type::numVertices; ++f )
                        pointIdInElt1[f] = elt1.point( f ).id();

                    index_type e0_0 = pointIdInElt0[elt0.e2p(idInElement0,0)];
                    index_type e0_1 = pointIdInElt0[elt0.e2p(idInElement0,1)];
                    index_type e1_0 = pointIdInElt1[elt1.e2p(idInElement1,0)];
                    index_type e1_1 = pointIdInElt1[elt1.e2p(idInElement1,1)];

                    if ( eltUseConnection0 )
                        continue;
                    if ( e0_0 != e1_0 )
                        elt.setEdgePermutation( idInElement0, edge_permutation_type( edge_permutation_type::REVERSE_PERMUTATION ) );
                }
            }
        }

    /**
     * update the entities of co-dimension 2
     */
    void updateEntitiesCoDimensionTwo() override
    {
        // no-op
    }

  private:
    friend class boost::serialization::access;
    template <class Archive>
    FEELPP_NO_EXPORT void serialize( Archive& ar, const unsigned int version )
    {
        ar& boost::serialization::base_object<super>( *this );
        DVLOG( 2 ) << "Serializing points\n";
        std::string str;
        str = "points";
        ar& str;
        ar& boost::serialization::base_object<super_points>( *this );
#if 1
        str = "faces";
        ar& str;
        DVLOG( 2 ) << "Serializing faces\n";
        ar& boost::serialization::base_object<super_faces>( *this );
#endif
        str = "elements";
        ar& str;
        DVLOG( 2 ) << "Serializing elements\n";
        ar& boost::serialization::base_object<super_elements>( *this );
    }
};

} // Feel

#endif /* __Mesh2D_H */
