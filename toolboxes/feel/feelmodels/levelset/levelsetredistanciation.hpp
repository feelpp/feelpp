//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
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
//! @file levelsetredistanciation.hpp
//! @author Thibaut Metivet <thibaut.metivet@gmail.com>
//! @date 25 Jul 2019
//! @copyright 2019 Feel++ Consortium
//!

#ifndef _LEVELSET_REDISTANCIATION_HPP
#define _LEVELSET_REDISTANCIATION_HPP 1

#include <memory>
#include <string>

#include <feel/feelvf/vf.hpp>

namespace Feel {

template<
    typename FunctionSpaceType
    >
class LevelSetRedistanciation
{
    public:
        //--------------------------------------------------------------------//
        // Typedefs
        typedef LevelSetRedistanciation<FunctionSpaceType> self_type;
        typedef std::shared_ptr<self_type> self_ptrtype;
        // Space
        typedef FunctionSpaceType functionspace_type;
        typedef std::shared_ptr<functionspace_type> functionspace_ptrtype;
        typedef typename functionspace_type::value_type value_type;
        // Element
        typedef typename functionspace_type::element_type element_type;
        typedef std::shared_ptr<element_type> element_ptrtype;
        // Mesh
        typedef typename functionspace_type::convex_type convex_type;
        static inline const uint16_type nDim = convex_type::nDim;
        static inline const uint16_type nOrderGeo = convex_type::nOrder;
        static inline const uint16_type nRealDim = convex_type::nRealDim;
        typedef typename functionspace_type::mesh_type mesh_type;
        typedef std::shared_ptr<mesh_type> mesh_ptrtype;
        // Periodicity
        typedef typename functionspace_type::periodicity_0_type periodicity_type;
        // Ranges
        typedef Range<mesh_type,MESH_ELEMENTS> range_elements_type;

        //--------------------------------------------------------------------//
        // Constructor/Destructor
        LevelSetRedistanciation( functionspace_ptrtype space, std::string const& prefix ) :
            M_space( space ),
            M_periodicity( boost::fusion::at_c<0>(space->periodicity()) ),
            M_prefix( prefix )
        {}

        // Prefix
        std::string const& prefix() const { return M_prefix; }

        // Space
        functionspace_ptrtype const& functionSpace() const { return M_space; }
        // Mesh
        mesh_ptrtype const& mesh() const { return this->functionSpace()->mesh(); }

        // Redistanciation
        virtual element_type run( element_type const& phi ) const =0;
        element_type operator()( element_type const& phi ) const { return this->run( phi ); }

    protected:
        std::string M_prefix;

        functionspace_ptrtype M_space;
        periodicity_type M_periodicity;
};

} // namespace Feel

#endif //_LEVELSET_REDISTANCIATION_HPP
