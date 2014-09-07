/*
This file is part of the Feel library
Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

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
/*! \file bareItems.h
\brief Special structures for handling mesh faces and sides
\version 0.0 Experimental   19/8/99. Luca Formaggia

Classes BareFace and BareEdge have been created to give an UNIQUE internal
representation for mesh faces and edges, allowing thus the construction of
Dof objects (which are naturally linked to mesh entities).

\par Introduction One of the paradigms chosen for the development of this
library is the fact that degrees of freedom (Dof) are linked to geometrical
entities.  Now if we have degrees of freedom associated, for instance, to
an Edge (like in a P2 Tetra) in order to build the global numbering of the
Dof and the association between local (element-wise) and global numbering,
I need to identify edges and give them a unique uint16_type. Yet, I may not want
want to build a full Edge (GeoElement2D) object: after all that I need is
the uint16_type of the edge and a way of computing the uint16_type's of the degrees of
freedom on the edge, all the remaining data of a full Edge object is not
necessarily needed.

Another related problem is how to uniquely identify a face or an edge in the mesh.


The dilemma has been resolved by creating the concept of a BareEdge and
BareFace (bare geometry items).  A bare geometry item is formed by the
minimal information required to uniquely identify it, namely 2
<tt>Point</tt>'s uint16_type 's for an edge and 3 <tt>Point</tt>'s uint16_type 's for the
Faces (it is enough also for Quad faces!). We build the bare items by
looping through the elements and obviously we make sure that the BareItem
uint16_type is consistent with that of the corresponding ``full item'' if the latter
has been instantiated.

Another <em>very important</em> issue is that of orientation.  There are
different ways of considerin orientation of a Face or an Edge. The first is
the <em>local</em> orientation of a Face or Egde of the reference finite
element. This is conventionally chosen when designing the finite
elements. For the faces, we have adopted the convention that the local
positive orientation is such that the face normal calculated with the right
hand rule is <em>outwardly</em> oriented. As for the edges, the local
orientation for a 3D element is more arbitrary.

However, for a 2D element, the positive orientation of an Edge is the one
which is in accordance with the right hand rule applied to that element.


When a Face or an edge is <em>active</em>, i.e. is effectively stored in
the mesh, then there is another obvious orientation, of global rather than
local nature: that induced by the way the edge or face is stored. For
boundary elements (faced in 3D or edges in 2D) it is compulsory that the
orientation of the stored item be consistent with the convention chosen for
the orientation of the domain boundary. More precisely, boundary elements
are stored so that the normal (always calculated following the righ hand
rule) is outward with respect to the domain.

However, there is the need of defining a <em>global</em> orientation also
for <em>non active</em> entities. This because we do not always store all
faces and all edges. We need then to choose a unique way to identify the
orientation of an Edge or of a Face <em>independently </em> from the fact
that they are active or not. We will call this orientation the <em>natural </em>
orientation. We have chosen the following
convention for natural orientation of faces and edges

<ul>
<li>The positive natural orientation of an  <em>Edge</em> is given by \f$V_{min} \rightarrow V_{max} \f$,
    \f$V_{min}\f$ being the Edge Vertex with smallest uint16_type</li>

<li>The positive natural orientation of a  <em>Face</em> is given by the cicle
    \f$V_{min} \rightarrow V_2\rightarrow V_3 \f$,
    \f$V_{min}\f$ being the Face Vertex with smallest uint16_type, \f$V_2\f$ the second smallest
    and \f$V_2\f$ the thirsd smallest.</li>
</ul>

Note that the latter definition applies both to triangular and to quad faces.

\warning If I want to associate boundary conditions I need the active
entity, since BareEdges do not store Marker data. (This is why the
RegionMesh classes treat boundary items in a rather special way).



\par Usage

Since BareItems are used only internally, we are (temporarily) omitting a
detailed documentation.
*/

#ifndef _MESHBAREITEMS_HH_
#define _MESHBAREITEMS_HH_
/*!
  Special routines to read meshes and special structures for
  sides and faces handling

  Classes BareFace and BareEdge have been created to give an UNIQUE
  Representation for mesh faces and edges and thus allow the construction
  of global tables or Fields.

*/
#include<utility>
#include<vector>
#include<map>
#include<algorithm>
#include<iostream>

#include <feel/feelcore/feel.hpp>

namespace Feel
{
/**
 * \class BarePoint
 *  \brief The Point basis class
 *  It contains the attributes common to all Points In particular, it
 * contains the one size_type's (first) of the point
 */
struct BarePoint
{
    /** Standard Constructor*/
    BarePoint() : first( 0 )
    {}

    /** constructor taking the size_type's */
    BarePoint( size_type i ) : first( i )
    {}

    /** First size_type which defines the Edge */
    size_type first;

};

/**
 * \class BareEdge
 *  \brief The Edge basis class
 *  It contains the attributes common to all Edges In particular, it
 * contains the two size_type's (first and second) of the points at the two ends
 * of the edge.
 * \invariant first < second
 */
struct BareEdge
{
    /** Standard Constructor*/
    BareEdge() : first( 0 ), second( 0 )
    {}

    /** constructor taking the size_type's */
    BareEdge( size_type i, size_type j ) : first( i ), second( j )
    {}

    /** First size_type which defines the Edge */
    size_type first;

    /** Second size_type which defines the Edge */
    size_type second;
};

/**
 * \class BareFace
 *  \brief The base Face class
 *
 *  Contains the common parameters of all Edges A Face contains three
 * size_type's (first, second and third) of points at face vertex, chosen so
 * to uniquely identify the face.
 * \invariant first<second<third.
 *
 */
struct BareFace
{
    /** Default constructor */
    BareFace() : first( 0 ), second( 0 ), third( 0 )
    {}

    /** Constructor that takes the size_type's as parameter */
    BareFace( size_type i, size_type j, size_type k ) : first( i ), second( j ), third( k )
    {}

    /**
     * Constructor that takes a BareEdge Object and an size_type. The
     * face is then identified by the size_type of the Points on the
     * BareEdge + <tt>i</tt>
     */
    BareFace( size_type i, const BareEdge & e ) : first( i ), second( e.first ), third( e.second )
    {}

    /** First size_type which defines the BareFace */
    size_type first;
    /** Second size_type which defines the BareFace */
    size_type second;
    /** Third size_type which defines the BareFace */
    size_type third;
};


    /** \defgroup BareItemsBuilder Global functions to build Bare Items.
        \ingroup Obsolet_Groups */
inline
std::pair<BarePoint, bool>
makeBarePoint( size_type const i )
{
    return std::make_pair( BarePoint( i ), true );
}
/**
 *  \ingroup BareItemsBuilder
 *
 *  \brief It creates a BareEdge end returns the orientation of the
 *  created edge with respect to the given data.
 *
 * \return is false if orientation  has been changed.
 * \param i is a Point size_type
 * \param j is a Point size_type
 *
 * The BareEdge that will be built is the one passing by <tt>i</tt>
 * and <tt>j</tt>. The orientation is i->j if the returned parameter
 * is a true.
 *
 * \pre i and j >0, i!=j
 */
inline
std::pair<BareEdge, bool>
makeBareEdge( size_type const i, size_type const j )
{
    if ( i < j )
    {
        return std::make_pair( BareEdge( i, j ), true );
    }

    else
    {
        return std::make_pair( BareEdge( j, i ), false );
        ;
    }
}
inline
std::pair<BareEdge, bool>
makeBareItem( size_type i, size_type j )
{
    return makeBareEdge( i, j );
}

/**
 * \ingroup BareItemsBuilder
 *
 * \param i is a Point size_type
 * \param j is a Point size_type
 * \brief It creates a BareEdge, ignoring orientation.
 *
 * The BareEdge that will be built is the one passing by <tt>i</tt> and <tt>j</tt>
 * A lighter version of MakeBareEdge, to be used
 * if orientation flag is not needed;
 * \pre i and j >0, i!=j
 */
inline
BareEdge
setBareEdge( size_type const i, size_type const j )
{
    BareEdge be;
    be.first = i < j ? i : j;
    be.second = ( i - be.first ) + j;
    return be;
}

/**
 * \ingroup BareItemsBuilder
 *
 * \param i is a Point size_type
 * \param j is a Point size_type
 *
 * \pre i and j >0
 * \brief It creates a non-standard BareEdge.
 *
 * Yet another lighter version of MakeBareEdge, without orientation,
 * To be used for non-oriented graphs.
 *
 * \warning It produces a BareEdge which does not comply with the
 * invariant of the class (first < second). It must be used only if
 * the BareEdge class is NOT used to uniquely identify edges.
 */
inline
BareEdge
setBareEdgeNo( size_type const i, size_type const j )
{
    return BareEdge( i, j );
}


/*! \ingroup BareItemsBuilder
 \brief It creates Bare Face objects from three Point size_type's
  \param  bool is false if orientation  has been changed.
  \param i is a Point size_type
  \param j is a Point size_type
  \param k is a Point size_type

  To be used for triangular faces.
  \pre i, j and k >0. i!=j!=k

*/
std::pair<BareFace, bool> makeBareFace( size_type i, size_type j, size_type k );
inline
std::pair<BareFace, bool>
makeBareItem( size_type i, size_type j, size_type k )
{
    return makeBareFace( i, j, k );
}

/*! \ingroup BareItemsBuilder
 \brief It creates Bare Face objects from four Point size_type's
  \param  bool is false if orientation  has been changed.
  \param i is a Point size_type
  \param j is a Point size_type
  \param k is a Point size_type
  \param l is a Point size_type

  To be used for triangular faces.
  \pre i, j and k >0. i!=j!=k

*/
std::pair<BareFace, bool> makeBareFace( size_type i, size_type j, size_type k, size_type l );

inline
std::pair<BareFace, bool>
makeBareItem( size_type i, size_type j, size_type k, size_type l )
{
    return makeBareFace( i, j, k, l );
}

/*! \defgroup comparison Comparison Operators
  \ingroup Obsolet_Groups
  Operators for comparing BareItems
*/

/*! \ingroup comparison
    inequality
*/
inline
bool
operator!=( const BareEdge & p1 , const BareEdge & p2 )
{
    return p1.first != p2.first || p1.second != p2.second;
}

/*! \ingroup comparison
     equality
*/
inline
bool
operator==( const BareEdge & p1 , const BareEdge & p2 )
{
    return p1.first == p2.first && p1.second == p2.second;
}

/*! \ingroup comparison
    greater than
*/
inline
bool
operator>( const BareEdge & e1 , const BareEdge & e2 )
{
    return e2.first > e1.first || ( e2.first == e1.first && e2.second > e1.second );
}

/*! \ingroup comparison
     greater-equal than
*/
inline
bool
operator>=( const BareEdge & e1 , const BareEdge & e2 )
{
    return e1 == e2 || e1 > e2;
}

/*! \ingroup comparison
    less than
*/
inline
bool
operator<( const BareEdge & e1 , const BareEdge & e2 )
{
    return e2.first < e1.first || ( e2.first == e1.first && e2.second < e1.second );
}

/*! \ingroup comparison
     less-equal than
*/
inline
bool
operator<=( const BareEdge & e1 , const BareEdge & e2 )
{
    return e1 == e2 || e1 < e2;
    ;
}

/*! \ingroup comparison*/
inline
bool
operator!=( const BareFace & p1 , const BareFace & p2 )
{
    return p1.first != p2.first || p1.second != p2.second || p1.third != p2.third;
}


/*! \ingroup comparison*/
inline
bool
operator==( const BareFace & p1 , const BareFace & p2 )
{
    return p1.first == p2.first && p1.second == p2.second && p1.third == p2.third;
}

/*! \ingroup comparison
  General functor for lexicographic comparison
*/
template <typename T>
struct cmpBareItem;

/*! \ingroup comparison
   Specialised functor for Edges
*/
// Specialisations
template <>
struct cmpBareItem<BarePoint> //!< The actual comparison operator
{
    bool operator() ( const BarePoint & e1, const BarePoint & e2 ) const
    {
        return e2.first > e1.first;
    }
};

// Specialisations
template <>
struct cmpBareItem<BareEdge> //!< The actual comparison operator
{
    bool operator() ( const BareEdge & e1, const BareEdge & e2 ) const
    {
        return e2.first > e1.first || ( e2.first == e1.first && e2.second > e1.second );
    }
};

/*! \ingroup comparison
  Specialised functor for Faces
*/
template <>
struct cmpBareItem<BareFace>
{
    bool operator() ( const BareFace & e1, const BareFace & e2 ) const
    {
        if ( e2.first > e1.first )
            return true;

        if ( e2.first == e1.first )
        {
            if ( e2.second > e1.second )
                return true;

            if ( e2.second == e1.second )
                return e2.third > e1.third;
        }

        return false;
    }
};

/**
 * \class BareItemsHandler
 * \brief Bare Items Handler
 *
 * This class handles mesh bare edges and faces construction. Used
 * only in mesh builders A BareItemsHandler is a specialisation of a
 * STL map which holds the pair formed by a bareitem and its
 * size_type.  The size_type is automatically generated if one uses
 * the method addIfNotThere
*/
template <typename BareItem>
class BareItemsHandler: public std::map<BareItem, Feel::size_type, cmpBareItem<BareItem> >
{
public:
    typedef std::map<BareItem, Feel::size_type, cmpBareItem<BareItem> > container;
    typedef typename container::size_type size_type;
    typedef typename container::iterator iterator;
    typedef typename container::const_iterator const_iterator;
    typedef std::pair<const BareItem, size_type> value_type;

    BareItemsHandler();

    //!< is the item there? I just ask
    bool isThere( BareItem const & ) const;

    //!< Returns size_type of a BareItem. 0 if not there
    size_type id( BareItem const & ) const;

    //!< To modify size_type of bareitem item in the list
    bool setId( BareItem const & item, size_type const i );

    //!< if not there adds it, the item size_type is autogenerated
    std::pair<size_type, bool> addIfNotThere( BareItem const & );

    //!<if not there adds it, and sets size_type id
    std::pair<size_type, bool> addIfNotThere( BareItem const &, const size_type id );

    //!< if it is there take it out (Id is lost)
    bool isThereDel( BareItem const & );

    //!< The # of entities ones actually stored.
    size_type howMany() const
    {
        return container::size();
    }

    //!< Max size_type currently in use
    size_type maxId() const
    {
        return M_id_count;
    }

    //!< Writes info in output
    void showMe() const;

private:
    size_type  M_id_count;
};

/*********************************************************************************
               IMPLEMENTATIONS
*********************************************************************************/
//
/*! \defgroup Helper Some helper functions
  \ingroup Obsolet_Groups */

//!\ingroup Helper
template <typename BareItem>
inline
size_type getId( std::pair<BareItem, size_type> const & i )
{
    return i.second;
}

//!\ingroup Helper
template <typename BareItem>
inline
BareItem getItem( std::pair<BareItem, size_type> const & i )
{
    return i.first;
}


//                   BareItemsHandler
template <class BareItem>
BareItemsHandler<BareItem>::BareItemsHandler()
    :
    M_id_count( 0 )
{ }


template <class BareItem>
inline
bool
BareItemsHandler<BareItem>::isThere( const BareItem & s ) const
{
    return find( s ) != container::end();
}

template <class BareItem>
inline
bool
BareItemsHandler<BareItem>::setId( const BareItem & s, size_type const id )
{
    const_iterator i = find( s );

    if ( i != container::end() )
    {
        i->second = id;
        return true;
    }

    else
    {
        return false;
    }

}

template <class BareItem>
inline
typename BareItemsHandler<BareItem>::size_type
BareItemsHandler<BareItem>::id( const BareItem & s ) const
{
    const_iterator i = find( s );

    if ( i != container::end() )
        return i->second;

    else
        return 0;
}

template <class BareItem>
inline
std::pair<typename BareItemsHandler<BareItem>::size_type, bool>
BareItemsHandler<BareItem>::addIfNotThere( const BareItem & s )
{
    std::pair<typename BareItemsHandler<BareItem>::iterator, bool> i( insert( std::make_pair( s, M_id_count ) ) );

    if ( i.second )
        ++M_id_count;

    return std::make_pair( ( i.first ) ->second, i.second );
}

template <class BareItem>
inline
std::pair<typename BareItemsHandler<BareItem>::size_type, bool>
BareItemsHandler<BareItem>::addIfNotThere( const BareItem & s, const size_type id )
{
    std::pair<typename BareItemsHandler<BareItem>::iterator, bool> i( insert( std::make_pair( s, id ) ) );
    // Set new id in any case.
    ( i.first ) ->second = id;

    // id count should grow +1
    //FEELPP_ASSERT( M_id_count == id )( M_id_count )( id ).error ( "invalid item id" );
    if ( i.second )
    {
        M_id_count = id;
        ++M_id_count;
    }

    // for consistency with other version.
    return std::make_pair( id, i.second );
}

template <class BareItem>
bool
BareItemsHandler<BareItem>::isThereDel( BareItem const & s )
{
    return erase( s ) != 0;
}

template <typename BareItem>
inline
void BareItemsHandler<BareItem>::showMe() const
{
    std::cout << "BareItemsHandler: " << std::endl;
    std::cout << "Number of Items stored: " << this->size() << std::endl;
    std::cout << "Max Id stored         : " << this->maxId() << std::endl;
    std::cout << "End of Information";
}
}
#endif
