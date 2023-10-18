/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2011-08-24

  Copyright (C) 2011 Université Joseph Fourier (Grenoble I)

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
   \file mesh0d.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2011-08-24
 */
#ifndef __Mesh0d_H
#define __Mesh0d_H 1


#include <iomanip>
#include <fstream>
#include <cstdlib>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/io.hpp>



#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/visitor.hpp>

#include <feel/feelmesh/meshbase.hpp>

#include <feel/feelmesh/geoelement.hpp>

#include <feel/feelmesh/elements.hpp>
#include <feel/feelmesh/points.hpp>

namespace Feel
{
/**
 *
 * \brief 0D mesh class
 *
 * \code
 * // create a 0D mesh made of simplex of order 1
 * Mesh0d<Simplex<0,1> > mesh;
 * \endcode
 *
 *  @author Christophe Prud'homme
 *  @see
 */
template<typename Shape, typename T = double, typename IndexT = uint32_type>
class Mesh0D
    :
        public VisitableBase<>,
        public MeshBase<IndexT>,
        public Elements<Shape,T>,
        public Points<Shape::nRealDim,T>
{
    // check at compilation time that the shape has indeed dimension 1
    BOOST_STATIC_ASSERT( Shape::nDim == 0 && Shape::nRealDim >= 1 );

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
    
    typedef Elements<Shape,T> super_elements;
    typedef typename super_elements::elements_type elements_type;
    typedef typename super_elements::element_type element_type;
    typedef typename super_elements::element_iterator element_iterator;
    typedef typename super_elements::element_const_iterator element_const_iterator;
    typedef typename super_elements::update_element_neighbor_type update_element_neighbor_type;
    using interprocess_element_iterator = element_iterator;
    using interprocess_element_const_iterator = element_const_iterator;
    using location_element_iterator = typename super_elements::location_element_iterator;
    using location_element_const_iterator = typename super_elements::location_element_const_iterator;
    using marker_element_iterator = typename super_elements::marker_element_iterator;
    using marker_element_const_iterator = typename super_elements::marker_element_const_iterator;
    using marker2_element_iterator = typename super_elements::marker2_element_iterator;
    using marker2_element_const_iterator = typename super_elements::marker2_element_const_iterator;
    using marker3_element_iterator = typename super_elements::marker3_element_iterator;
    using marker3_element_const_iterator = typename super_elements::marker3_element_const_iterator;

    typedef super_elements super_faces;
    typedef elements_type faces_type;

    typedef Points<nRealDim,T> super_points;
    typedef typename super_points::points_type points_type;
    typedef typename super_points::point_type point_type;

    typedef Mesh0D<Shape,T> self_type;
    typedef std::shared_ptr<self_type> self_ptrtype;

    using face_type = point_type;
    using face_iterator = element_iterator;
    using face_const_iterator = element_const_iterator;
    using edge_type = point_type;
    using edge_iterator = element_iterator;
    using edge_const_iterator = element_const_iterator;

    using pid_face_iterator = face_iterator;
    using pid_face_const_iterator = face_const_iterator;
    using location_face_iterator = location_element_iterator;
    using location_face_const_iterator = location_element_const_iterator;
    using interprocess_face_iterator = location_element_iterator;
    using interprocess_face_const_iterator = location_element_const_iterator;
    using marker_face_iterator = marker_element_iterator;
    using marker_face_const_iterator = marker_element_const_iterator;
    using marker2_face_iterator = marker2_element_iterator;
    using marker2_face_const_iterator = marker2_element_const_iterator;
    using marker3_face_iterator = marker3_element_iterator;
    using marker3_face_const_iterator = marker3_element_const_iterator;

    using pid_edge_iterator = face_iterator;
    using pid_edge_const_iterator = face_const_iterator;
    using location_edge_iterator = location_element_iterator;
    using location_edge_const_iterator = location_element_const_iterator;
    using interprocess_edge_iterator = location_element_iterator;
    using interprocess_edge_const_iterator = location_element_const_iterator;
    using marker_edge_iterator = marker_element_iterator;
    using marker_edge_const_iterator = marker_element_const_iterator;
    using marker2_edge_iterator = marker2_element_iterator;
    using marker2_edge_const_iterator = marker2_element_const_iterator;
    using marker3_edge_iterator = marker3_element_iterator;
    using marker3_edge_const_iterator = marker3_element_const_iterator;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    /**
     * default constructor
     */
    explicit Mesh0D( std::string const& name = "", worldcomm_ptr_t const& worldComm = Environment::worldCommPtr() )
        :
        super_visitable(),
        super( name, 0, nRealDim, worldComm ),
        super_elements( worldComm ),
        super_points( worldComm )
        {}

    Mesh0D( Mesh0D const& m ) = default;
    Mesh0D( Mesh0D && m ) = default;

    ~Mesh0D() override
        {}

    //@}

    /** @name Operator overloads
     */
    //@{

    Mesh0D& operator=( Mesh0D const& m ) = default;
    Mesh0D& operator=( Mesh0D && m ) = default;

    //@}

    /** @name Accessors
     */
    //@{

    /**
     * \return \p true if all containers are empty, \p false otherwise
     */
    bool isEmpty() const
        {
            return ( super_elements::isEmpty() &&
                     super_points::isEmpty() );
        }


    /**
     * \return the number of elements
     */
    size_type numElements() const
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
    size_type numFaces() const
        {
            return 0;
        }


    /**
     * \return the number of points
     */
    size_type numPoints() const
        {
            return this->points().size();
        }

    //@}

    /** @name  Mutators
     */
    //@{

    void setWorldComm( worldcomm_ptr_t const& _worldComm ) override
        {
            MeshBase<IndexT>::setWorldComm( _worldComm );
            this->setWorldCommPoints( _worldComm );
        }

    //@}

    /** @name  Methods
     */
    //@{


    /**
     * clear out all data from the mesh, \p isEmpty() should return
     * \p true after a \p clear()
     */
    virtual void clear()
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
    void renumber()
        {
            FEELPP_ASSERT( 0 ).error( "invalid call" );
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
    void updateEntitiesCoDimensionTwo()
        {
            // no-op
        }
private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
        {
            ar & boost::serialization::base_object<super_elements>( *this );
            ar & boost::serialization::base_object<super_points>( *this );
        }



};




} // Feel

#endif /* __Mesh0D_H */
