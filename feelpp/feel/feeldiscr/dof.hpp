/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 2010-04-20

 Copyright (C) 2013-2016 Feel++ Consortium

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

#include <tuple>

namespace Feel
{
/**
 * @brief class that represents a degree of freedom
 *
 * @see DofTable, FaceDof
 */
template<typename SizeT=uint32_type>
class Dof 
    :
        public std::tuple<SizeT>
        //public std::tuple<size_type, int16_type, bool>
//std::tuple<size_type, int16_type, bool, uint16_type, bool, size_type>
{
    //typedef std::tuple<size_type, int16_type, bool, uint16_type, bool, size_type> super;
    //typedef std::tuple<size_type, int16_type, bool> super;
    typedef std::tuple<SizeT> super;
public:


    /** @name Constants
     */
    //@{


    //@}

    /** @name Typedefs
     */
    //@{

    using size_type = SizeT;
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
            std::get<0>(*this) =  gid;
            //this->get<1>() =  1;
            //this->get<2>() =  false;

        }
    Dof( size_type gid, int16_type s )
        :
        super( )
        {
            std::get<0>(*this) =  gid;
            //this->get<1>() =  s;
            //this->get<2>() =  false;

        }

    Dof( std::tuple<size_type, int16_type, bool> const& t )
        :
        super( )
        {
            std::get<0>(*this) =  std::get<0>(t);
            //this->get<1>() =  t.get<1>();
            //this->get<2>() =  t.get<2>();

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
        super(_index ) //, _sign, per)
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
            std::get<0>(*this) =  t;
            //this->get<1>() =  1;
            //this->get<2>() =  false;
            return *this;
        }
    Dof& operator=( std::tuple<size_type, int16_type, bool> const& t )
        {
            std::get<0>(*this) =  std::get<0>(t);
            //this->get<1>() =  t.get<1>();
            //this->get<2>() =  t.get<2>();
            return *this;
        }

    //@}

    /** @name Accessors
     */
    //@{

    /// @return the global index
    size_type index() const
        {
            return std::get<0>(*this);
        }

    /// @return the sign
    int16_type sign() const
        {
            return 1;//this->get<1>();
        }
    /// @return if periodic
    bool isPeriodic() const
        {
            return false;//this->get<2>();
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
            std::get<0>(*this) = id;
        }

    /**
     * set the global dof
     */
    void set( size_type _index, int16_type _sign, bool per )
        {
            std::get<0>(*this) =  _index;
            //this->get<1>() =  _sign;
            //this->get<2>() =  per; 
        }


    //@}


    /** @name  Methods
     */
    //@{


    //@}
};

template<typename SizeT=uint32_type>
inline
std::ostream&
operator<<( std::ostream& __os, Dof<SizeT> const& __dof )
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
template<typename SizeT=uint32_type>
struct FaceDof : public std::tuple<SizeT, int16_type, bool, uint16_type, uint16_type>
{
    typedef std::tuple<SizeT, int16_type, bool, uint16_type, uint16_type> super;
public:

    using size_type = SizeT;
    
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
            std::get<0>(*this) =  gid;
            std::get<1>(*this) =  1;
            std::get<2>(*this) =  false;
            std::get<3>(*this) =  -1;


        }

    FaceDof( std::tuple<size_type, int16_type, bool> const& t )
        :
        super( )
        {
            std::get<0>(*this) =  std::get<0>(t);
            std::get<1>(*this) =  std::get<1>(t);
            std::get<2>(*this) =  std::get<2>(t);

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
        super(_index, _sign, per, ld, -1 )
        {
        }
    FaceDof( Dof<size_type> const& d, uint16_type ldinface, uint16_type ldinelt   )
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
            std::get<0>(*this) =  t;
            std::get<1>(*this) =  1;
            std::get<2>(*this) =  false;
            return *this;
        }
    FaceDof& operator=( std::tuple<size_type, int16_type, bool> const& t )
        {
            std::get<0>(*this) =  std::get<0>(t);
            std::get<1>(*this) =  std::get<1>(t);
            std::get<2>(*this) =  std::get<2>(t);
            return *this;
        }

    //@}

    /** @name Accessors
     */
    //@{

    /// @return the global index
    size_type index() const
        {
            return std::get<0>(*this);
        }

    /// @return the sign
    int16_type sign() const
        {
            return std::get<1>(*this);
        }
    /// @return if periodic
    bool isPeriodic() const
        {
            return std::get<2>(*this);
        }

    /// @return the local dof in the element
    uint16_type localDof() const
        {
            return std::get<3>(*this);
        }
    /// @return the local dof in the face
    uint16_type localDofInFace() const
        {
            return std::get<4>(*this);
        }
    /// @return the local dof in the face
    uint16_type localDofInEntity() const
        {
            return std::get<4>(*this);
        }
    /// @return the local dof in element
    uint16_type localDofInElement() const                                \
        {
            return std::get<3>(*this);
        }

    //@}

    /** @name  Mutators
     */
    //@{
    // set the global dof id
    void setIndex( size_type id )
        {
            std::get<0>(*this) = id;
        }

    //@}


    /** @name  Methods
     */
    //@{


    //@}

};
//@ Alias for FaceDof
template<typename SizeT=uint32_type>
using EntityDof =  FaceDof<SizeT>;

template<int NC = 1, typename SizeT=uint32_type>
class LocalDof: public std::pair<SizeT,uint16_type>
{
public:
    using size_type = SizeT;
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

template<int NC,typename SizeT=uint32_type>
std::ostream&
operator<<( std::ostream& __os, LocalDof<NC,SizeT> const& __dof )
{
    __os << "-----------Dof-Info------------\n"
         << "elementId             : " << __dof.elementId() << "\n"
         << "localDof              : " << __dof.localDof() << "\n"
         << "localDofPerComponent  : " << __dof.localDofPerComponent() << "\n"
         << "nComponents : " << __dof.nComponents() << "\n";
    return __os;
}

template<int NC = 1,typename SizeT = uint32_type>
class LocalDofSet : public std::vector<LocalDof<NC,SizeT>>
{
  public:
    using size_type = SizeT;
    typedef std::vector<LocalDof<NC,size_type>> super;
    typedef LocalDof<NC,size_type> localdof_type;
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
    LocalDofSet const& update( size_type eid, uint16_type nLocalDof )
    {
        if ( nLocalDof == this->size() )
            return this->update( eid );
        this->resize( nLocalDof );
        for(uint16_type i = 0; i < nLocalDof; ++i )
            this->at( i ) = localdof_type( eid, i );
        return *this;
    }

};

} // Feel
#endif /* FEELPP_DOF_HPP */
