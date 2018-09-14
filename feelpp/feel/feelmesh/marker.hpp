/*
 This file is part of the Feel library
 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano
 Copyright (C) 2006 Université Joseph Fourier (UJF)

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
/*!
 * \file markers.hpp
 * \brief A Simple implementation of Markers
 * \author Christophe Prud'homme
 *
 * this is based on Luca Formaggia marker class in Feel
 */
#ifndef HH_MARKERS_HH_
#define HH_MARKERS_HH_

#include <limits>

#include <boost/serialization/serialization.hpp>

#include <feel/feelcore/feel.hpp>
#include <boost/detail/identifier.hpp>

namespace Feel
{
/// \cond detail
class Marker1 : public boost::detail::identifier< size_type, Marker1 >
{
public:
    typedef boost::detail::identifier< size_type, Marker1 >::value_type value_type;
    Marker1()                           : boost::detail::identifier<size_type,Marker1>( 0 ) {}
    explicit Marker1( value_type v )    : boost::detail::identifier<size_type,Marker1>( v ) {}
    Marker1 & operator=( value_type v )
    {
        this->assign( v );
        return *this;
    }
    bool isOn() const { return value() != 0; }
    bool isOff() const { return value() == 0; }

private:

    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void save( Archive & ar, const unsigned int version ) const
        {
            size_type v = this->value();
            ar & v;
        }

    template<class Archive>
    void load( Archive & ar, const unsigned int version )
        {
            size_type v;
            ar &  v;
            this->assign( v );
        }

    BOOST_SERIALIZATION_SPLIT_MEMBER()

};

class Marker2 : public boost::detail::identifier< size_type, Marker2 >
{
public:
    typedef boost::detail::identifier< size_type, Marker2 >::value_type value_type;
    Marker2()                           : boost::detail::identifier<size_type,Marker2>( 0 ) {}
    explicit Marker2( value_type v )    : boost::detail::identifier<size_type,Marker2>( v ) {}
    Marker2 & operator=( value_type v )
    {
        this->assign( v );
        return *this;
    }
    bool isOn() const { return value() != 0; }
    bool isOff() const { return value() == 0; }

private:

    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void save( Archive & ar, const unsigned int version ) const
        {
            size_type v = this->value();
            ar & v;
        }

    template<class Archive>
    void load( Archive & ar, const unsigned int version )
        {
            size_type v;
            ar &  v;
            this->assign( v );
        }

    BOOST_SERIALIZATION_SPLIT_MEMBER()

};
class Marker3 : public boost::detail::identifier< size_type, Marker3 >
{
public:
    typedef boost::detail::identifier< size_type, Marker3 >::value_type value_type;
    Marker3()                           : boost::detail::identifier<size_type,Marker3>( 0 ) {}
    explicit Marker3( value_type v )    : boost::detail::identifier<size_type,Marker3>( v ) {}
    Marker3 & operator=( value_type v )
    {
        this->assign( v );
        return *this;
    }
    bool isOn() const { return value() != 0; }
    bool isOff() const { return value() == 0; }

private:

    friend class boost::serialization::access;
    // When the class Archive corresponds to an output archive, the
    // & operator is defined similar to <<.  Likewise, when the class Archive
    // is a type of input archive the & operator is defined similar to >>.
    template<class Archive>
    void save( Archive & ar, const unsigned int version ) const
        {
            size_type v = this->value();
            ar & v;
        }

    template<class Archive>
    void load( Archive & ar, const unsigned int version )
        {
            size_type v;
            ar &  v;
            this->assign( v );
        }

    BOOST_SERIALIZATION_SPLIT_MEMBER()

};


namespace detail
{
struct by_marker {};
struct by_marker2 {};
struct by_marker3 {};
struct by_interprocessdomain {};
struct by_location {};
struct by_pid {};
struct by_element {};
struct by_entity {};
struct by_ghostcell {};
}
/// \endcond detail
} // Feel
#endif

