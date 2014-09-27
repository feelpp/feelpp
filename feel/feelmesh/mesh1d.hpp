/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-11-09

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2007-2010 Université Joseph Fourier (Grenoble I)

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
   \file mesh1d.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-11-09
 */
#ifndef __Mesh1D_H
#define __Mesh1D_H 1


#include <iomanip>
#include <fstream>
#include <cstdlib>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/multi_array.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/mem_fun.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/numeric/ublas/io.hpp>



#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/visitor.hpp>

#include <feel/feelmesh/meshbase.hpp>

#include <feel/feelmesh/geoelement.hpp>

#include <feel/feelmesh/elements.hpp>
#include <feel/feelmesh/faces.hpp>
#include <feel/feelmesh/points.hpp>

namespace Feel
{
/**
 * \class Mesh1D
 * \brief 1D mesh class
 *
 * \code
 * // create a 1D mesh made of simplex of order 1
 * Mesh1D<Simplex<1,1> > mesh;
 *
 * // create a 1D mesh made of simplex of order 2
 * Mesh1D<Simplex<1,2> > mesh;
 * \endcode
 *
 *  @author Christophe Prud'homme
 *  @see
 */
template<typename Shape>
class Mesh1D
    :
        public VisitableBase<>,
        public MeshBase,
        public Elements<Shape>,
        public Points<Shape::nRealDim>,
        public Faces<typename Shape::template shape<0,Shape::nOrder,Shape::nRealDim>::type,
                     typename Elements<Shape>::element_type>
{
    // check at compilation time that the shape has indeed dimension 1
    BOOST_STATIC_ASSERT( Shape::nDim == 1 );

public:


    /** @name Typedefs
     */
    //@{

    static const uint16_type nDim = Shape::nRealDim;
    static const uint16_type nRealDim = Shape::nRealDim;

    typedef typename VisitableBase<>::return_type return_type;

    typedef VisitableBase<> super_visitable;
    typedef MeshBase super;

    typedef Elements<Shape> super_elements;
    typedef typename super_elements::elements_type elements_type;
    typedef typename super_elements::element_type element_type;
    typedef typename super_elements::element_iterator element_iterator;
    typedef typename super_elements::element_const_iterator element_const_iterator;
    typedef typename super_elements::update_element_neighbor_type update_element_neighbor_type;

    typedef Points<nRealDim> super_points;
    typedef typename super_points::points_type points_type;
    typedef typename super_points::point_type point_type;

    typedef Faces<typename Shape::template shape<0,Shape::nOrder,Shape::nRealDim>::type,typename super_elements::element_type> super_faces;
    typedef typename super_faces::face_iterator face_iterator;
    typedef typename super_faces::face_const_iterator face_const_iterator;
    typedef typename super_faces::faces_type faces_type;
    typedef typename super_faces::face_type face_type;
    typedef face_type edge_type;

    typedef super_faces super_edges;
    typedef typename super_edges::marker_face_iterator marker_edge_iterator;
    typedef typename super_edges::marker_face_const_iterator marker_edge_const_iterator;
    typedef typename super_edges::location_face_iterator location_edge_iterator;
    typedef typename super_edges::location_face_const_iterator location_edge_const_iterator;
    typedef typename super_edges::pid_face_iterator pid_edge_iterator;
    typedef typename super_edges::pid_face_const_iterator pid_edge_const_iterator;

    typedef Mesh1D<Shape> self_type;
    typedef boost::shared_ptr<self_type> self_ptrtype;

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
    Mesh1D( WorldComm const& worldComm = Environment::worldComm() )
        :
        super_visitable(),
        super( worldComm ),
        super_elements( worldComm ),
        super_points( worldComm ),
        super_faces( worldComm )
        {}

/**
 * copy constructor
 */
    Mesh1D( Mesh1D const & m )
        :
        super_visitable(),
        super( m ),
        super_elements( m ),
        super_points( m ),
        super_faces( m )
        {}

/**
 * destructor
 */
    ~Mesh1D()
        {}

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
                super_faces::operator=( m );
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
    bool isEmpty() const
        {
            return ( super_elements::isEmpty() &&
                     super_points::isEmpty() &&
                     super_faces::isEmpty() );
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
            return this->faces().size();
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


//@}

/** @name  Methods
 */
//@{

    virtual void setWorldComm( WorldComm const& _worldComm )
        {
            this->setWorldCommMeshBase( _worldComm );
            this->setWorldCommElements( _worldComm );
            this->setWorldCommFaces( _worldComm );
            this->setWorldCommPoints( _worldComm );
        }

/**
 * clear out all data from the mesh, \p isEmpty() should return
 * \p true after a \p clear()
 */
    virtual void clear()
        {
            this->elements().clear();
            this->points().clear();
            this->faces().clear();
            FEELPP_ASSERT( isEmpty() ).error( "all mesh containers should be empty after a clear." );
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
            ar & boost::serialization::base_object<super>( *this );
            DVLOG(2) << "Serializing points\n";
            ar & boost::serialization::base_object<super_points>( *this );
            DVLOG(2) << "Serializing faces\n";
            ar & boost::serialization::base_object<super_faces>( *this );
            DVLOG(2) << "Serializing elements\n";
            ar & boost::serialization::base_object<super_elements>( *this );
        }


};




} // Feel

#endif /* __Mesh1D_H */
