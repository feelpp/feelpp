/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 2010-04-20

 Copyright (C) 2013 Feel++ Consortium

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
 \file dof.hpp
 \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 \date 2010-04-20
 */
#ifndef FEELPP_DOF_HPP
#define FEELPP_DOF_HPP 1

#include <boost/tuple/tuple.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/composite_key.hpp>
#include <boost/multi_index/mem_fun.hpp>
#include <boost/multi_index/ordered_index.hpp>

namespace Feel
{
/**
 * @brief class that represents a degree of freedom
 *
 * @see DofTable, FaceDof
 */
class Dof 
    : 
        public boost::tuple<size_type, int16_type, bool>
//boost::tuple<size_type, int16_type, bool, uint16_type, bool, size_type>
{
    //typedef boost::tuple<size_type, int16_type, bool, uint16_type, bool, size_type> super;
    typedef boost::tuple<size_type, int16_type, bool> super;
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
    Dof()
        :
        super()
        {
        }

    Dof( size_type gid )
        :
        super( )
        {
            this->get<0>() =  gid;
            this->get<1>() =  1;
            this->get<2>() =  false;

        }
    Dof( size_type gid, int16_type s )
        :
        super( )
        {
            this->get<0>() =  gid;
            this->get<1>() =  s;
            this->get<2>() =  false;

        }

    Dof( boost::tuple<size_type, int16_type, bool> const& t )
        :
        super( )
        {
            this->get<0>() =  t.get<0>();
            this->get<1>() =  t.get<1>();
            this->get<2>() =  t.get<2>();

        }


    /**
     *
     */
#if 0
    Dof( size_type _index, int16_type _sign, bool per, uint16_type _entity = 0, bool _location = false, size_type _marker = 0  )
        :
        super(_index, _sign, per, _entity, _location, _marker)
        {
        }
#else
    Dof( size_type _index, int16_type _sign, bool per )
        :
        super(_index, _sign, per)
        {
        }
#endif
    //! copy constructor
    Dof( Dof const & dof )
        :
        super( dof )
        {}

    /// move constructor
    Dof( Dof && d ) = default;

    //! destructor
    ~Dof()
        {}

    //@}

    /** @name Operator overloads
     */
    //@{

    //! copy operator
    Dof& operator=( Dof const & o ) = default;

    /// move operator=
    Dof& operator=( Dof&& d ) = default;

    Dof& operator=( size_type t )
        {
            this->get<0>() =  t;
            this->get<1>() =  1;
            this->get<2>() =  false;
            return *this;
        }
    Dof& operator=( boost::tuple<size_type, int16_type, bool> const& t )
        {
            this->get<0>() =  t.get<0>();
            this->get<1>() =  t.get<1>();
            this->get<2>() =  t.get<2>();
            return *this;
        }

    //@}

    /** @name Accessors
     */
    //@{

    /// @return the global index
    size_type index() const
        {
            return this->get<0>();
        }

    /// @return the sign
    int16_type sign() const
        {
            return this->get<1>();
        }
    /// @return if periodic
    bool isPeriodic() const
        {
            return this->get<2>();
        }
#if 0
    /// @return the entity type (0: vertex, 1:edge, 2:face, 3:volume)
    uint16_type entity() const
        {
            return this->get<3>();
        }

    /// @return the location
    bool isOnBoundary() const
        {
            return this->get<4>();
        }

    /// @return the marker
    size_type marker() const
        {
            return this->get<5>();
        }
#endif

    //@}

    /** @name  Mutators
     */
    //@{
    // set the global dof id
    void setIndex( size_type id )
        {
            this->get<0>() = id;
        }

    /**
     * set the global dof
     */
    void set( size_type _index, int16_type _sign, bool per )
        {
            this->get<0>() =  _index;
            this->get<1>() =  _sign;
            this->get<2>() =  per; 
        }


    //@}


    /** @name  Methods
     */
    //@{


    //@}
};


inline
std::ostream&
operator<<( std::ostream& __os, Dof const& __dof )
{
    __os << "-----------Dof-Info------------\n"
         << "index        : " << __dof.index() << "\n"
         << "sign         : " << __dof.sign() << "\n"
         << "isPeriodic   : " << __dof.isPeriodic() << "\n";
#if 0
    << "isOnBoundary : " << __dof.isOnBoundary() << "\n"
    << "marker       : " << __dof.marker() << "\n";
#endif
    return __os;
}

/**
 * @brief Describe a Dof on a Face
 * 
 * @details
 * the data structure is a tuple containing
 *  - global dof id
 *  - sign of the dof
 *  - is dof periodic
 *  - local dof id in element
 *  - local dof id in face
 */
struct FaceDof : public boost::tuple<size_type, int16_type, bool, uint16_type, uint16_type>
{
    typedef boost::tuple<size_type, int16_type, bool, uint16_type, uint16_type> super;
public:

    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    FaceDof()
        :
        super()
        {
        }

    FaceDof( size_type gid )
        :
        super( )
        {
            this->get<0>() =  gid;
            this->get<1>() =  1;
            this->get<2>() =  false;
            this->get<3>() =  -1;


        }

    FaceDof( boost::tuple<size_type, int16_type, bool> const& t )
        :
        super( )
        {
            this->get<0>() =  t.get<0>();
            this->get<1>() =  t.get<1>();
            this->get<2>() =  t.get<2>();

        }
    FaceDof( super const& t )
        :
        super( t )
        {
        }


    /**
     *
     */
    FaceDof( size_type _index, int16_type _sign, bool per, uint16_type ld  )
        :
        super(_index, _sign, per, ld )
        {
        }
    FaceDof( Dof const& d, uint16_type ldinface, uint16_type ldinelt   )
        :
        super(d.index(), d.sign(), d.isPeriodic(), ldinelt, ldinface )
        {
        }

    //! copy constructor
    FaceDof( FaceDof const & dof )
        :
        super( dof )
        {}

    //! destructor
    ~FaceDof()
        {}

    //@}

    /** @name Operator overloads
     */
    //@{

    //! copy operator
    FaceDof& operator=( FaceDof const & o )
        {
            if ( this != &o )
            {
                super::operator=( o );
            }

            return *this;
        }    

    FaceDof& operator=( FaceDof && o ) = default;

    FaceDof& operator=( size_type t )
        {
            this->get<0>() =  t;
            this->get<1>() =  1;
            this->get<2>() =  false;
            return *this;
        }
    FaceDof& operator=( boost::tuple<size_type, int16_type, bool> const& t )
        {
            this->get<0>() =  t.get<0>();
            this->get<1>() =  t.get<1>();
            this->get<2>() =  t.get<2>();
            return *this;
        }

    //@}

    /** @name Accessors
     */
    //@{

    /// @return the global index
    size_type index() const
        {
            return this->get<0>();
        }

    /// @return the sign
    int16_type sign() const
        {
            return this->get<1>();
        }
    /// @return if periodic
    bool isPeriodic() const
        {
            return this->get<2>();
        }

    /// @return the local dof in the element
    uint16_type localDof() const
        {
            return this->get<3>();
        }
    /// @return the local dof in the face
    uint16_type localDofInFace() const
        {
            return this->get<4>();
        }
    /// @return the local dof in element
    uint16_type localDofInElement() const                                \
        {
            return this->get<3>();
        }

    //@}

    /** @name  Mutators
     */
    //@{
    // set the global dof id
    void setIndex( size_type id )
        {
            this->get<0>() = id;
        }

    //@}


    /** @name  Methods
     */
    //@{


    //@}

};
#if 0
typedef multi_index::multi_index_container<
    Dof,
    multi_index::indexed_by<

        // sort by less<int> on index()
        multi_index::ordered_unique<multi_index::const_mem_fun<Dof,
                                                               size_type,
                                                               &Dof::index> >,

        // sort by less<int> on entity()
        multi_index::ordered_non_unique<multi_index::tag<detail::by_entity>,
                                        multi_index::const_mem_fun<Dof,
                                                                   uint16_type,
                                                                   &Dof::entity> >,
        // sort by less<int> on location()
        multi_index::ordered_non_unique<multi_index::tag<detail::by_location>,
                                        multi_index::const_mem_fun<Dof,
                                                                   bool,
                                                                   &Dof::isOnBoundary> >,
        // sort by less<int> on marker()
        multi_index::ordered_non_unique<multi_index::tag<detail::by_marker>,
                                        multi_index::const_mem_fun<Dof,
                                                                   Marker1,
                                                                   &Dof::marker> >
        > > dof_container_type;
#endif

template<int NC = 1>
class LocalDof: public std::pair<size_type,uint16_type>
{
public:
    typedef std::pair<size_type,uint16_type> super;

    static constexpr uint16_type nComponents() { return NC; }

    LocalDof()
        :
        super( std::make_pair( 0, 0 ) )
        {}
    LocalDof( size_type e )
        :
        super( std::make_pair( e, 0 ) )
        {}
    LocalDof( size_type e, uint16_type l )
        :
        super( std::make_pair( e, l ) )
        {}
    LocalDof( std::pair<int,int> const& p )
        :
        super( p )
        {}
    size_type elementId() const { return this->first; }
    uint16_type localDof() const { return this->second; }
    uint16_type localDofPerComponent() const { return this->second/nComponents(); }
    // returns the local dof component given the number of local dof per component @arg nLocalDofPerComponent
    uint16_type component( uint16_type nLocalDofPerComponent ) const { return this->second/nLocalDofPerComponent; }
    
    void setLocalDof( uint16_type l ) { this->second = l; }
    void set( size_type e, uint16_type l ) { this->first=e; this->second = l; }

};

template<int NC>
std::ostream&
operator<<( std::ostream& __os, LocalDof<NC> const& __dof )
{
    __os << "-----------Dof-Info------------\n"
         << "elementId             : " << __dof.elementId() << "\n"
         << "localDof              : " << __dof.localDof() << "\n"
         << "localDofPerComponent  : " << __dof.localDofPerComponent() << "\n"
         << "nComponents : " << __dof.nComponents() << "\n";
    return __os;
}

template<int NC = 1>
class LocalDofSet : public std::vector<LocalDof<NC>>
{
  public:
    typedef std::vector<LocalDof<NC>> super;
    typedef LocalDof<NC> localdof_type;
    static constexpr uint16_type nComponents() { return NC; }
    LocalDofSet()
        :
        super()
    {}

    LocalDofSet( size_type eid, uint16_type nLocalDof )
        :
        super(nLocalDof)
    {
        for(uint16_type i = 0; i < nLocalDof; ++i )
        {
            this->at( i ) = localdof_type( eid, i );
        }
    }
    LocalDofSet const& update( size_type eid )
    {
        DCHECK( !this->empty() ) << "Invalid Local Dof Set";
        std::for_each( this->begin(), this->end(), [eid]( localdof_type& d ) { d.first = eid; } );
        return *this;
    }
};

} // Feel
#endif /* FEELPP_DOF_HPP */
