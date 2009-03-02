/*
 This file is part of the Life library
 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano
 Copyright (C) 2006 Université Joseph Fourier (UJF)

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
/*!
 * \file markers.hpp
 * \brief A Simple implementation of Markers
 * \author Christophe Prud'homme
 *
 * this is based on Luca Formaggia marker class in Life
 */
#ifndef HH_MARKERS_HH_
#define HH_MARKERS_HH_

#include <limits>

#include <life/lifecore/life.hpp>
#include <boost/detail/identifier.hpp>

namespace Life
{
#if 0
/**
   @class Marker
   @brief Base marker class.

   It stores an integral type, aliased to flag_type, which may be
   used for marking a geometrical entity. The typical use is to
   specify boundary conditions or material properties associated with
   the entity.  The actual boundary conditions will be handled in the
   Dof class. During the creation of a field, the markers are
   post-processed to furnish the correct boundary conditions.
*/
class Marker
{
public:
    //! Null flag as default
    explicit Marker()
        :
        flag( invalid_flag_type_value )
    {}

    //! Give the flag
    explicit Marker( flag_type & p )
        :
        flag( p )
    {}

    // copy constructor
    Marker( Marker const & m )
        :
        flag( m.flag )
    {}

    //! Extract marker flag
    flag_type marker() const { return flag; }

    //! Returns the null flag (cannot be modified)
    flag_type const & invalidFlag() const { return invalid_flag_type_value; }

    //! Set marker to given value
    flag_type set( flag_type const & c ) { flag = c; return flag; }

    //! Set marker to given value only if unset
    flag_type update( flag_type const & c );

    //! Sets the flag to the stronger flag of two given markers
    flag_type setStronger( flag_type const & p1, flag_type const & p2 );

    //! Sets the flag to the weaker flag of two given markers
    flag_type setWeaker( flag_type const & p1, flag_type const & p2 );

    //! If marker flag is unset, is sets it to that of the argument, otherwise
    //! is sets it to  the stronger flag between the stored one
    //! and the one provided by the argument.
    flag_type setStronger( flag_type const & p );

    //! If marker flag is unset , it sets it to that of the argument, otherwise
    //! is sets it to  the weaker flag between the stored one
    //! and the one provided by the argument.
    flag_type setWeaker( flag_type const & p );

    //! It enquires if marker flag is different than the nullflag
    bool isSet() const { return flag != invalidFlag(); }

    //! It enquires if marker flag is different than the nullflag
    bool isUnset() const { return flag == invalidFlag(); }

    //! Put marker to nullflag
    void unset() { flag = invalidFlag(); }

    //! Helper function that prints a marker Flag
    std::ostream & printFlag( flag_type const f, std::ostream & out ) const;

    //! Helper function that prints "this" marker flag
    std::ostream & printFlag( std::ostream & out ) const;

        /**
     * @brief Selects the stronger between two flags
     *
     * A dimensional geometric entity G_i may inherit the stronger flag
     * among adiacent geometric entities of greater dimensions.  For
     * example a boundary point Point with an unset Flag may inherit the
     * strongerst Flag of the adjacent boundary faces.
     * It returns invalid_flag_type_value if any of the entity a or b is a invalid_flag_type_value.
     */
    static flag_type strongerFlag( flag_type const & a, flag_type const & b )
    {
        return a > b ? a : b ;
    }

    /**
     * @brief Selects the weaker between two flags
     *
     *  A lower dimensional geometric entity G_i may inherit the weaker
     *  flag among the Vertices of the entity.  For example a boundary
     *  Face with an unset Flag may inherit the weakest Flag of its
     *  Vertices. This method may also be used to attribute the Flag to
     *  Nodes generated on high order elements.
     *  It returns invalid_flag_type_value if any of the entity a or b is a invalid_flag_type_value.
     */
    static flag_type weakerFlag( flag_type const & a, flag_type const & b )
    {
        if(a==invalid_flag_type_value)return b;
        if(b==invalid_flag_type_value)return a;
        return a < b ? a : b ;
    }

protected:
    flag_type flag;
};

inline
flag_type
Marker::update( flag_type const & c )
{
    if ( flag == invalidFlag() )
        return set( c );
}

inline
flag_type
Marker::setStronger( flag_type const & p1, flag_type const & p2 )
{
    return set( Marker::strongerFlag( p1, p2 ) );
}

inline
flag_type
Marker::setWeaker( flag_type const & p1, flag_type const & p2 )
{
    return set( Marker::weakerFlag( p1, p2 ) );
}

inline
flag_type
Marker::setStronger( flag_type const & p )
{
    if ( isUnset() )
        return flag = p;
    return set( Marker::strongerFlag( this->marker(), p ) );
}

inline
flag_type
Marker::setWeaker( flag_type const & p )
{
    if ( isUnset() )
        return flag = p;
    return set( Marker::weakerFlag( this->marker(), p ) );
}

inline
std::ostream &
Marker::printFlag( flag_type const f, std::ostream & out ) const
{
    if ( f == invalidFlag() )
        out << "UNSET";
    else
        out << f;
    return out;
}
inline
std::ostream &
Marker::printFlag( std::ostream & out ) const
{
    return printFlag( flag, out );
}
#else



class Marker1 : public boost::detail::identifier< int64_type, Marker1 >
{
public:
    typedef boost::detail::identifier< int64_type, Marker1 >::value_type value_type;
    Marker1()                           : boost::detail::identifier<int64_type,Marker1>(0){}
    explicit Marker1( value_type v )    : boost::detail::identifier<int64_type,Marker1>(v){}
    Marker1 & operator=( value_type v ) { this->assign(v); return *this; }
};

class Marker2 : public boost::detail::identifier< int64_type, Marker2 >
{
public:
    typedef boost::detail::identifier< int64_type, Marker2 >::value_type value_type;
    Marker2()                           : boost::detail::identifier<int64_type,Marker2>(0){}
    explicit Marker2( value_type v )    : boost::detail::identifier<int64_type,Marker2>(v){}
    Marker2 & operator=( value_type v ) { this->assign(v); return *this; }
};
class Marker3 : public boost::detail::identifier< int64_type, Marker3 >
{
public:
    typedef boost::detail::identifier< int64_type, Marker3 >::value_type value_type;
    Marker3()                           : boost::detail::identifier<int64_type,Marker3>(0){}
    explicit Marker3( value_type v )    : boost::detail::identifier<int64_type,Marker3>(v){}
    Marker3 & operator=( value_type v ) { this->assign(v); return *this; }
};
#endif
namespace detail
{
struct by_marker{};
struct by_marker2{};
struct by_marker3{};
struct by_interprocessdomain{};
struct by_location{};
struct by_pid{};
struct by_element{};
}
} // Life
#endif

