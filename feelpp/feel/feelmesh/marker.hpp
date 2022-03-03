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
#ifndef FEELPP_MESH_MARKERS_H
#define FEELPP_MESH_MARKERS_H

//#include <limits>
#include <boost/serialization/serialization.hpp>
#include <feel/feelcore/feel.hpp>
//#include <boost/detail/identifier.hpp>

namespace Feel
{

template <typename ValueType = flag_type>
class Marker : public std::set<ValueType>
{
    using super_type = std::set<ValueType>;
public:
    using value_type = ValueType;

    Marker() = default;
    explicit Marker( value_type v ) : super_type({ v }) {}
    Marker( super_type const& s ) : super_type(s) {}
    Marker( Marker const& ) = default;
    Marker( Marker && ) = default;

    Marker & operator=( Marker const& ) = default;
    Marker & operator=( Marker && ) = default;

    bool isOn() const { return !this->empty(); }
    bool isOff() const { return this->empty(); }

    void assign( value_type v ) { super_type::operator=({v}); }

    value_type value() const { CHECK( !this->empty() ) << "no marker"; return *this->begin(); }

private:

    friend class boost::serialization::access;

    template <class Archive>
    void serialize( Archive& ar, const unsigned int version )
    {
        ar& boost::serialization::base_object<super_type>( *this );
    }
};


} // Feel
#endif

