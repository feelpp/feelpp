//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! Author(s) : 
//!     Thibaut Metivet <thibaut.metivet@inria.fr>
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
//! @file fastmarchingdofstatus.hpp
//! @author Thibaut Metivet <thibaut.metivet@inria.fr>
//! @date 12 Dec 2019
//! @copyright 2019 Feel++ Consortium
//! @copyright 2019 INRIA
//!

#ifndef _FASTMARCHINGDOFSTATUS_HPP
#define _FASTMARCHINGDOFSTATUS_HPP 1

namespace Feel {

struct FastMarchingDofStatus
{
    typedef uint8_type status_type;
    
    static constexpr status_type FAR = 1 << 0;
    static constexpr status_type CLOSE = 1 << 1;
    static constexpr status_type DONE = 1 << 2;

    static constexpr status_type FIX = 1 << 3;
    static constexpr status_type OLD = 1 << 4;
    static constexpr status_type NEW = 1 << 5;

    static constexpr status_type CLOSE_OLD = CLOSE | OLD;
    static constexpr status_type CLOSE_NEW = CLOSE | NEW;

    static constexpr status_type DONE_FIX = DONE | FIX;
    static constexpr status_type DONE_OLD = DONE | OLD;
    static constexpr status_type DONE_NEW = DONE | NEW;

    FastMarchingDofStatus() = default;
    FastMarchingDofStatus( status_type const& s ): status(s) {}
    operator status_type() const { return status; }

    bool operator&( status_type const s ) const
    {
        return ( this->status & s );
    }
    FastMarchingDofStatus & operator&=( status_type const s )
    {
        this->status &= s;
        return *this;
    }

    bool operator|( status_type const s ) const
    {
        return ( this->status | s );
    }
    FastMarchingDofStatus & operator|=( status_type const s )
    {
        this->status |= s;
        return *this;
    }

    bool operator^( status_type const s ) const
    {
        return ( this->status ^ s );
    }

    FastMarchingDofStatus operator~() const
    {
        return FastMarchingDofStatus( ~( this->status ) );
    }

    status_type status;
};

}

namespace boost {
namespace serialization {
    template<class Archive>
    void serialize( Archive & ar, Feel::FastMarchingDofStatus & s, const unsigned int version )
    {
        ar & s.status;
    }
}
namespace mpi {
    template <>
    struct is_mpi_datatype<Feel::FastMarchingDofStatus> : mpl::true_ { };
}
}

#endif

