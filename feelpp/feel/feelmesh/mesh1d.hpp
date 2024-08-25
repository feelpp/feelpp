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
//! This file provides the header of the 1D mesh data structure
//!
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 05 Feb 2017
//! @copyright 2005,2006 EPFL
//! @copyright 2007-2010 Université Joseph Fourier (Grenoble I)
//! @copyright 2011-2017 Feel++ Consortium
//!
#ifndef FEELPP_MESH1D_HPP
#define FEELPP_MESH1D_HPP 1

#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <fstream>
#include <iomanip>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>


#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/shared_ptr.hpp>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/visitor.hpp>

#include <feel/feelmesh/meshbase.hpp>

#include <feel/feelmesh/geoelement.hpp>

#include <feel/feelmesh/elements.hpp>
#include <feel/feelmesh/faces.hpp>
#include <feel/feelmesh/points.hpp>

namespace Feel
{
//!
//! @brief 1D mesh class
//! @ingroup Mesh
//!
//! @code
//! // create a 1D mesh made of simplex of order 1
//! Mesh1D<Simplex<1,1> > mesh;
//!
//! // create a 1D mesh made of simplex of order 2
//! Mesh1D<Simplex<1,2> > mesh;
//! @endcode
//!
//! @author Christophe Prud'homme
//! @see Mesh2D, Mesh3D
//!
template<typename Shape, typename T = double, typename IndexT = uint32_type>
class Mesh1D
    :
        public VisitableBase<>,
        public MeshBase<IndexT>,
        public Elements<Shape,T,IndexT>,
        public Points<Shape::nRealDim,T,IndexT, SubFaceOf<typename Elements<Shape,T,IndexT>::element_type> >
{
    // check at compilation time that the shape has indeed dimension 1
    BOOST_STATIC_ASSERT( Shape::nDim == 1 );

  public:
    /** @name Typedefs
     */
    //@{

    static inline const uint16_type nDim = Shape::nRealDim;
    static inline const uint16_type nRealDim = Shape::nRealDim;

    typedef typename VisitableBase<>::return_type return_type;

    typedef VisitableBase<> super_visitable;
    typedef MeshBase<IndexT> super;

    using index_type = typename super::index_type;
    using size_type = typename super::size_type;

    typedef Elements<Shape,T,IndexT> super_elements;
    typedef typename super_elements::elements_type elements_type;
    typedef typename super_elements::element_type element_type;
    typedef typename super_elements::element_iterator element_iterator;
    typedef typename super_elements::element_const_iterator element_const_iterator;
    typedef typename super_elements::update_element_neighbor_type update_element_neighbor_type;

    typedef Points<nRealDim,T,IndexT,SubFaceOf<typename Elements<Shape,T,IndexT>::element_type>> super_points;
    typedef typename super_points::points_type points_type;
    typedef typename super_points::point_type point_type;

    typedef super_points super_faces;
    typedef typename super_faces::point_iterator face_iterator;
    typedef typename super_faces::point_const_iterator face_const_iterator;
    typedef typename super_faces::points_type faces_type;
    typedef typename super_faces::point_type face_type;
    typedef typename super_faces::points_reference_wrapper_type faces_reference_wrapper_type;
    typedef typename super_faces::points_reference_wrapper_ptrtype faces_reference_wrapper_ptrtype;
    typedef typename faces_reference_wrapper_type::iterator face_reference_wrapper_iterator;
    typedef typename faces_reference_wrapper_type::const_iterator face_reference_wrapper_const_iterator;

    typedef face_type edge_type;

    typedef super_faces super_edges;
    typedef typename super_edges::point_iterator edge_iterator;
    typedef typename super_edges::point_const_iterator edge_const_iterator;
    typedef typename super_edges::points_reference_wrapper_type edges_reference_wrapper_type;
    typedef typename super_edges::points_reference_wrapper_ptrtype edges_reference_wrapper_ptrtype;
    typedef typename edges_reference_wrapper_type::iterator edge_reference_wrapper_iterator;
    typedef typename edges_reference_wrapper_type::const_iterator edge_reference_wrapper_const_iterator;

    typedef Mesh1D<Shape,T, IndexT> self_type;

    typedef std::shared_ptr<self_type> self_ptrtype;

    typedef typename element_type::vertex_permutation_type vertex_permutation_type;
    typedef typename element_type::edge_permutation_type edge_permutation_type;

    typedef typename super::face_processor_type face_processor_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    /**
     * default constructor
     */
    explicit Mesh1D( std::string const& name = "", worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() )
        : super_visitable(),
          super( name, 1, nRealDim, worldComm ),
          super_elements( worldComm ),
          super_points( worldComm )
    {
    }

    //! copy constructor
    Mesh1D( Mesh1D const& m )
        : super_visitable(),
          super( m ),
          super_elements( m ),
          super_points( m )
    {
    }

    //! destructor
    ~Mesh1D() override
    {
    }

    //@}

    /** @name Operator overloads
     */
    //@{

    Mesh1D& operator=( Mesh1D const& m )
    {
        if ( this != &m )
        {
            super::operator=( m );
            super_elements::operator=( m );
            super_points::operator=( m );
        }

        return *this;
    }

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
                 super_points::isEmpty() );
    }

    //!
    //! @brief get the number of elements in the mesh
    //! @return the number of elements in the mesh
    //!
    size_type numElements() const override
    {
        return this->elements().size();
    }

    /**
     * \return the number of faces in an element
     */
    size_type numLocalFaces() const
    {
        return super_elements::element_type::numLocalFaces;
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
        return this->points().size();
    }
    size_type numEdges() const
    {
        return this->points().size();
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

    faces_type & faces()
    {
        return this->points();
    }

    faces_type const& faces() const
    {
        return this->points();
    }

    bool hasFace( size_type i ) const
    {
        return this->hasPoint( i );
    }

    face_type const& face( size_type i ) const
    {
        return this->point( i );
    }

    face_const_iterator faceIterator( size_type i ) const
    {
        return this->pointIterator( i );
    }
    face_iterator faceIterator( size_type i )
    {
        return this->pointIterator( i );
    }
    face_const_iterator faceIterator( face_type const& face ) const
    {
        return faceIterator( face.id() );
    }
    face_iterator faceIterator( face_type const& face )
    {
        return faceIterator( face.id() );
    }

    face_iterator beginFace()
    {
        return this->beginPoint();
    }
    face_const_iterator beginFace() const
    {
        return this->beginPoint();
    }
    face_iterator endFace()
    {
        return this->endPoint();
    }
    face_const_iterator endFace() const
    {
        return this->endPoint();
    }

     std::pair<face_iterator,bool> addFace( face_type& f )
        {
            return this->addPoint( f );
        }
    std::pair<face_iterator,bool> addFace( face_type&& f )
        {
             return this->addPoint( f );
        }

     face_iterator eraseFace( face_iterator it )
        {
            return this->erasePoint( it );
        }

    std::tuple<face_reference_wrapper_const_iterator,face_reference_wrapper_const_iterator,faces_reference_wrapper_ptrtype>
    internalFaces( rank_type p = invalid_rank_type_value ) const
        {
            return this->internalPoints( p );
        }
    std::tuple<face_reference_wrapper_const_iterator,face_reference_wrapper_const_iterator,faces_reference_wrapper_ptrtype>
    facesOnBoundary( rank_type p = invalid_rank_type_value ) const
        {
            return this->boundaryPoints( p );
        }

    std::tuple<face_reference_wrapper_const_iterator,face_reference_wrapper_const_iterator,faces_reference_wrapper_ptrtype>
    facesWithMarkerByType( uint16_type markerType, std::set<flag_type> const& markerFlags, rank_type p = invalid_rank_type_value ) const
        {
            return this->pointsWithMarkerByType( markerType, markerFlags, p );
        }
     std::tuple<face_reference_wrapper_const_iterator,face_reference_wrapper_const_iterator,faces_reference_wrapper_ptrtype>
    facesWithMarkerByType( uint16_type markerType, flag_type m, rank_type p = invalid_rank_type_value ) const
        {
            return this->pointsWithMarkerByType( markerType, m, p );
        }
    std::tuple<face_reference_wrapper_const_iterator,face_reference_wrapper_const_iterator,faces_reference_wrapper_ptrtype>
    facesWithMarker( flag_type m = invalid_flag_type_value, rank_type p = invalid_rank_type_value ) const
        {
            return this->facesWithMarkerByType( 1, m, p );
        }
    std::tuple<face_reference_wrapper_const_iterator,face_reference_wrapper_const_iterator,faces_reference_wrapper_ptrtype>
    facesWithMarker2( flag_type m = invalid_flag_type_value, rank_type p = invalid_rank_type_value ) const
        {
            return this->facesWithMarkerByType( 2, m, p );
        }
    std::tuple<face_reference_wrapper_const_iterator,face_reference_wrapper_const_iterator,faces_reference_wrapper_ptrtype>
    facesWithMarker3( flag_type m = invalid_flag_type_value, rank_type p = invalid_rank_type_value ) const
        {
            return this->facesWithMarkerByType( 3, m, p );
        }

    std::tuple<face_reference_wrapper_const_iterator,face_reference_wrapper_const_iterator,faces_reference_wrapper_ptrtype>
    facesWithProcessId( rank_type p = invalid_rank_type_value ) const
        {
            return this->pointsWithProcessId( p );
        }
    std::tuple<face_reference_wrapper_const_iterator,face_reference_wrapper_const_iterator,faces_reference_wrapper_ptrtype>
    interProcessFaces( rank_type neighbor_pid = invalid_rank_type_value ) const
        {
            return this->interProcessPoints( neighbor_pid );
        }
    
    

    void setWorldComm( worldcomm_ptr_t const& _worldComm ) override
    {
        MeshBase<IndexT>::setWorldComm( _worldComm );
        this->setWorldCommElements( _worldComm );
        this->setWorldCommPoints( _worldComm );
    }

    /**
     * clear out all data from the mesh, \p isEmpty() should return
     * \p true after a \p clear()
     */
    void clear() override
        {
            super::clear();
            super_elements::clear();
            super_points::clear();
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
        CHECK( false ) << "invalid call";
    }

    /**
     * update permutation of entities of co-dimension 1
     */
    void updateEntitiesCoDimensionOnePermutation() 
    {
        // no-op
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

    template<class Archive>
    FEELPP_NO_EXPORT void serialize( Archive & ar, const unsigned int version )
        {
            ar & boost::serialization::base_object<super>( *this );
            DVLOG(2) << "Serializing points\n";
            ar & boost::serialization::base_object<super_points>( *this );
            DVLOG(2) << "Serializing elements\n";
            ar & boost::serialization::base_object<super_elements>( *this );
        }

};

} // Feel

#endif /* __Mesh1D_H */
