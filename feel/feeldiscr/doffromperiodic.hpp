/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <prudhomme@unistra.fr>
       Date: 2013-10-23

  Copyright (C) 2013 Université de Strasbourg

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
   \file doffromperiodic.hpp
   \author Christophe Prud'homme <prudhomme@unistra.fr>
   \date 2013-10-23
 */
#ifndef FEELPP_DofFromPeriodic_H
#define FEELPP_DofFromPeriodic_H 1

namespace Feel
{
/**
 * \class DofFromPeriodic
 * \brief brief description
 *
 * @author Christophe Prud'homme
 * @see
 */
class DofFromPeriodic
{
public:


  /** @name Constants
*/
 //@{


//@}

/** @name Typedefs
*/
 //@{


//@}

/** @name Constructors, destructor
*/
//@{

//! default constructor
DofFromPeriodic();
  //! copy constructor
DofFromPeriodic( DofFromPeriodic const & );
  //! destructor
~DofFromPeriodic();

//@}

/** @name Operator overloads
*/
//@{

//! copy operator
DofFromPeriodic& operator=( DofFromPeriodic const & o)
{
 if (this != &o )
 {
 }
 return *this;
}
//@}

/** @name Accessors
*/
//@{


//@}

/** @name  Mutators
*/
//@{


//@}

/** @name  Methods
*/
//@{


//@}



protected:

private:
    /**
     * Add a new periodic dof to the dof map and to the list of periodic dof \p
     * periodic_dof for a given tag \p tag
     *
     * \arg __elt geometric element
     * \arg __face  face of the element that holds the periodic dof
     * \arg next_free_dof integer that will get possibilty incremented and holds the current number of dof
     * \arg periodic_dof table of periodic dof
     * \arg tag the boundary tag associated with the periodic dof
     */
    void addVertexPeriodicDof( element_type const& __elt,
                               face_type const& __face,
                               size_type& next_free_dof,
                               std::map<size_type,periodic_dof_map_type>& periodic_dof,
                               size_type tag )
    {
        addVertexPeriodicDof( __elt, __face, next_free_dof, periodic_dof, tag, mpl::bool_<(fe_type::nDofPerVertex>0)>() );
    }
    void addVertexPeriodicDof( element_type const& __elt,face_type const& __face,size_type& next_free_dof,std::map<size_type,periodic_dof_map_type>& periodic_dof,size_type tag, mpl::bool_<false> ) {}
    void addVertexPeriodicDof( element_type const& __elt,face_type const& __face,size_type& next_free_dof,std::map<size_type,periodic_dof_map_type>& periodic_dof,size_type tag, mpl::bool_<true> );

    void addEdgePeriodicDof( element_type const& __elt,
                             face_type const& __face,
                             size_type& next_free_dof,
                             std::map<size_type,periodic_dof_map_type>& periodic_dof,
                             size_type tag )
    {
        static const bool cond = fe_type::nDofPerEdge > 0;
        addEdgePeriodicDof( __elt, __face, next_free_dof, periodic_dof, tag , mpl::bool_<cond>(), mpl::int_<nDim>() );
    }
    void addEdgePeriodicDof( element_type const& __elt,face_type const& __face,size_type& next_free_dof,std::map<size_type,periodic_dof_map_type>& periodic_dof, size_type tag, mpl::bool_<false>, mpl::int_<1> ) {}
    void addEdgePeriodicDof( element_type const& __elt,face_type const& __face,size_type& next_free_dof,std::map<size_type,periodic_dof_map_type>& periodic_dof, size_type tag, mpl::bool_<false>, mpl::int_<2> ) {}
    void addEdgePeriodicDof( element_type const& __elt,face_type const& __face,size_type& next_free_dof,std::map<size_type,periodic_dof_map_type>& periodic_dof, size_type tag, mpl::bool_<false>, mpl::int_<3> ) {}
    void addEdgePeriodicDof( element_type const& __elt,face_type const& __face,size_type& next_free_dof,std::map<size_type,periodic_dof_map_type>& periodic_dof, size_type tag, mpl::bool_<true>, mpl::int_<1> ) {}
    void addEdgePeriodicDof( element_type const& __elt,face_type const& __face,size_type& next_free_dof,std::map<size_type,periodic_dof_map_type>& periodic_dof, size_type tag, mpl::bool_<true>, mpl::int_<2> );
    void addEdgePeriodicDof( element_type const& __elt,face_type const& __face,size_type& next_free_dof,std::map<size_type,periodic_dof_map_type>& periodic_dof, size_type tag, mpl::bool_<true>, mpl::int_<3> );


    void addFacePeriodicDof( element_type const& __elt,face_type const& __face,size_type& next_free_dof,std::map<size_type,periodic_dof_map_type>& periodic_dof,size_type tag )
    {
        static const bool cond = fe_type::nDofPerFace > 0;
        addFacePeriodicDof( __elt, __face, next_free_dof, periodic_dof, tag, mpl::bool_<cond>() );
    }
    void addFacePeriodicDof( element_type const& __elt,face_type const& __face,size_type& next_free_dof,std::map<size_type,periodic_dof_map_type>& periodic_dof,size_type tag, mpl::bool_<false> ) {}
    void addFacePeriodicDof( element_type const& __elt,face_type const& __face,size_type& next_free_dof,std::map<size_type,periodic_dof_map_type>& periodic_dof,size_type tag, mpl::bool_<true> );

};

template<typename MeshType, typename FEType, typename PeriodicityType>
void
DofTable<MeshType, FEType, PeriodicityType>::addVertexPeriodicDof( element_type const& __elt,
        face_type const& __face,
        size_type& next_free_dof,
        std::map<size_type,periodic_dof_map_type>& periodic_dof,
        size_type tag,
        mpl::bool_<true> )
{
    // store the element and local dof id for further
    // reference when inserting the associated global dof
    // id of the element adjacent to the face

    size_type iElAd = __face.ad_first();
    FEELPP_ASSERT( iElAd != invalid_size_type_value )( __face.id() ).error( "[periodic]invalid face/element in face" );
    Feel::detail::ignore_unused_variable_warning( iElAd );

    // local id of the face in its adjacent element
    uint16_type iFaEl = __face.pos_first();
    FEELPP_ASSERT( iFaEl != invalid_uint16_type_value ).error ( "invalid element index in face" );

    // loop on face vertices
    for ( uint16_type iVeFa = 0; iVeFa < face_type::numVertices; ++iVeFa )
    {
        // local vertex number (in element)
        uint16_type iVeEl = element_type::fToP( iFaEl, iVeFa );
        Feel::detail::ignore_unused_variable_warning( iVeEl );

        FEELPP_ASSERT( iVeEl != invalid_uint16_type_value ).error( "invalid local dof" );

        // Loop number of Dof per vertex
        for ( uint16_type l = 0; l < fe_type::nDofPerVertex; ++l )
        {
            uint16_type lid = iVeEl * fe_type::nDofPerVertex + l;
            //const size_type gDof = global_shift + ( __elt.point( i ).id() ) * fe_type::nDofPerVertex + l;
            const size_type gDof = ( __elt.point( iVeEl ).id() ) * fe_type::nDofPerVertex + l;

            VLOG(2) << "add vertex periodic doc " << next_free_dof << " in element " << __elt.id() << " lid = " << lid << "\n";
            size_type dof_id = next_free_dof;
            // next_free_dof might be incremented if a new dof is created
            bool inserted = this->insertDof( __elt.id(), lid, iVeEl, boost::make_tuple( 0, 0, gDof ), 0, next_free_dof, 1, true, 0, __elt.point( iVeEl ).marker()  );
            VLOG(2) << "vertex periodic dof inserted : " << inserted << "\n";

            const int ncdof = is_product?nComponents:1;
            for ( int c1 = 0; c1 < ncdof; ++c1 )
            {
                // add the pair (elt, lid) to the map associated
                // with dof_id, one dof can be shared by several
                // elements
                dof_id = localToGlobal( __elt.id(), lid, c1 ).index();
                periodic_dof[tag].insert( std::make_pair( dof_id, boost::make_tuple( __elt.id(), lid, c1, gDof, 0 ) ) );

                VLOG(2) << "added vertex periodic dof " <<  __elt.id() << ", " <<  lid << ", " << boost::get<0>( localToGlobal( __elt.id(), lid, c1 ) ) << "\n";
            }
        }

    }

}

template<typename MeshType, typename FEType, typename PeriodicityType>
void
DofTable<MeshType, FEType, PeriodicityType>::addEdgePeriodicDof( element_type const& __elt,
        face_type const& __face,
        size_type& next_free_dof,
        std::map<size_type,periodic_dof_map_type>& periodic_dof,
        size_type tag,
        mpl::bool_<true>,
        mpl::int_<2> )
{
#if 0
    // id of the element adjacent to the face
    // \warning NEED TO INVESTIGATE THIS
    size_type iElAd = __face.ad_first();
    FEELPP_ASSERT( iElAd != invalid_size_type_value )
    ( __face.id() ).error( "[DofTable::buildBoundaryDof] invalid face/element in face" );
#endif // 0

    // local id of the face in its adjacent element
    uint16_type iFaEl = __face.pos_first();
    FEELPP_ASSERT( iFaEl != invalid_uint16_type_value ).error ( "invalid element index in face" );
#if !defined(NDEBUG)
    DVLOG(4) << " local face id : " << iFaEl << "\n";
#endif

    // Loop number of DofTable per edge
    for ( uint16_type l = 0; l < fe_type::nDofPerEdge; ++l )
    {
        uint16_type lid = element_type::numVertices*fe_type::nDofPerVertex + iFaEl * fe_type::nDofPerEdge + l;
        //const size_type gDof = global_shift + ( __elt.point( i ).id() ) * fe_type::nDofPerVertex + l;
        const size_type gDof = ( __elt.edge( iFaEl ).id() ) * fe_type::nDofPerEdge + l;

        DVLOG(4) << "add edge periodic dof " << next_free_dof << " in element " << __elt.id() << " lid = " << lid << "\n";
        size_type dof_id = next_free_dof;
        // next_free_dof might be incremented if a new dof is created
        bool inserted = this->insertDof( __elt.id(), lid, iFaEl, boost::make_tuple( 1, 0, gDof ), 0, next_free_dof, 1, true, 0, __elt.edge( iFaEl ).marker() );
        DVLOG(4) << "edge periodic dof inserted (1 or 0) : " << inserted << "\n";

        const int ncdof = is_product?nComponents:1;
        for ( int c1 = 0; c1 < ncdof; ++c1 )
        {
            // add the pair (elt, lid) to the map associated
            // with dof_id, one dof can be shared by several
            // elements
            dof_id = boost::get<0>( localToGlobal( __elt.id(), lid, c1 ) );
            periodic_dof[tag].insert( std::make_pair( dof_id, boost::make_tuple( __elt.id(), lid, c1, gDof, 1 ) ) );

            DVLOG(4) << "added edge periodic dof " <<  __elt.id() << ", " <<  lid << ", " << boost::get<0>( localToGlobal( __elt.id(), lid, c1 ) ) << "\n";
        }

    }
}
template<typename MeshType, typename FEType, typename PeriodicityType>
void
DofTable<MeshType, FEType, PeriodicityType>::addEdgePeriodicDof( element_type const& __elt,
        face_type const& __face,
        size_type& next_free_dof,
        std::map<size_type,periodic_dof_map_type>& periodic_dof,
        size_type tag,
        mpl::bool_<true>,
        mpl::int_<3> )
{
    //BOOST_STATIC_ASSERT( face_type::numEdges );

    // id of the element adjacent to the face
    // \warning NEED TO INVESTIGATE THIS
    size_type iElAd = __face.ad_first();
    Feel::detail::ignore_unused_variable_warning( iElAd );
    FEELPP_ASSERT( iElAd != invalid_size_type_value )( __face.id() ).error( "[Dof::buildBoundaryDof] invalid face/element in face" );

    // local id of the face in its adjacent element
    uint16_type iFaEl = __face.pos_first();
    FEELPP_ASSERT( iFaEl != invalid_uint16_type_value ).error ( "invalid element index in face" );
#if !defined(NDEBUG)
    DVLOG(4) << " local face id : " << iFaEl << "\n";
#endif

    // loop on face vertices
    for ( uint16_type iEdFa = 0; iEdFa < face_type::numEdges; ++iEdFa )
    {
        // local edge number (in element)
        uint16_type iEdEl = element_type::fToE( iFaEl, iEdFa );
        Feel::detail::ignore_unused_variable_warning( iEdEl );
        FEELPP_ASSERT( iEdEl != invalid_uint16_type_value ).error( "invalid local dof" );

        // Loop number of Dof per edge
        for ( uint16_type l = 0; l < fe_type::nDofPerEdge; ++l )
        {

            uint16_type lid = element_type::numVertices*fe_type::nDofPerVertex + iFaEl * fe_type::nDofPerEdge + l;
            //const size_type gDof = global_shift + ( __elt.point( i ).id() ) * fe_type::nDofPerVertex + l;
            size_type gDof = __elt.edge( l ).id() * fe_type::nDofPerEdge;

            if ( __elt.edgePermutation( l ).value()  == edge_permutation_type::IDENTITY )
            {
                gDof += l ; // both nodal and modal case
            }

            else if ( __elt.edgePermutation( l ).value()  == edge_permutation_type::REVERSE_PERMUTATION )
            {
                gDof += fe_type::nDofPerEdge - 1 - l ;
            }

            DVLOG(4) << "add periodic doc " << next_free_dof << " in element " << __elt.id() << " lid = " << lid << "\n";
            size_type dof_id = next_free_dof;
            // next_free_dof might be incremented if a new dof is created
            bool inserted = this->insertDof( __elt.id(), lid, l, boost::make_tuple( 1, 0, gDof ), 0, next_free_dof, 1, true, 0, __elt.edge( l ).marker() );
            DVLOG(4) << "periodic dof inserted : " << inserted << "\n";

            const int ncdof = is_product?nComponents:1;
            for ( int c1 = 0; c1 < ncdof; ++c1 )
            {
                // add the pair (elt, lid) to the map associated
                // with dof_id, one dof can be shared by several
                // elements
                dof_id = boost::get<0>( localToGlobal( __elt.id(), lid, c1 ) );
                periodic_dof[tag].insert( std::make_pair( dof_id, boost::make_tuple( __elt.id(),  lid, c1, gDof, 1 ) ) );

                DVLOG(4) << "added " <<  __elt.id() << ", " <<  lid << ", " << boost::get<0>( localToGlobal( __elt.id(), lid, c1 ) ) << "\n";
            }
        }
    }

}
template<typename MeshType, typename FEType, typename PeriodicityType>
void
DofTable<MeshType, FEType, PeriodicityType>::addFacePeriodicDof( element_type const& __elt,
        face_type const& __face,
        size_type& next_free_dof,
        std::map<size_type,periodic_dof_map_type>& periodic_dof,
        size_type tag,
        mpl::bool_<true> )
{
#if 0
    // id of the element adjacent to the face
    // \warning NEED TO INVESTIGATE THIS
    size_type iElAd = __face.ad_first();
    FEELPP_ASSERT( iElAd != invalid_size_type_value )( __face.id() ).error( "[Dof::buildBoundaryDof] invalid face/element in face" );

    // local id of the face in its adjacent element
    uint16_type iFaEl = __face.pos_first();
    FEELPP_ASSERT( iFaEl != invalid_uint16_type_value ).error ( "invalid element index in face" );
#if !defined(NDEBUG)
    DVLOG(2) << " local face id : " << iFaEl << "\n";
#endif

    // Loop on number of Dof per face
    for ( uint16_type l = 0; l < fe_type::nDofPerFace; ++l )
    {
        auto temp = this->localToGlobal( iElAd,
                                         element_type::numVertices*fe_type::nDofPerVertex +
                                         element_type::numEdges*fe_type::nDofPerEdge +
                                         iFaEl * fe_type::nDofPerFace + l,
                                         c );
        M_face_l2g[ __face_it->id()][ lc++ ] = boost::make_tuple( boost::get<0>( temp ),boost::get<1>( temp ),boost::get<2>( temp ),
                                                element_type::numVertices*fe_type::nDofPerVertex +
                                                element_type::numEdges*fe_type::nDofPerEdge +
                                                iFaEl * fe_type::nDofPerFace + l );
    }

#endif // 0
}

}
#endif /* FEELPP_DofFromPeriodic_H */
