/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
 
 This file is part of the Feel++ library
 
 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 13 Sep 2015
 
 Copyright (C) 2015 Feel++ Consortium
 
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
#ifndef FEELPP_POLY_DOF_HPP
#define FEELPP_POLY_DOF_HPP 1

#include <vector>
#include <feel/feelcore/feel.hpp>

namespace Feel {

class FiniteElementDof 
{
public:
    using super = std::vector<uint16_type>;
    FiniteElementDof() = default;
    FiniteElementDof(FiniteElementDof const&) = default;
    FiniteElementDof(FiniteElementDof &&) = default;
    
    constexpr FiniteElementDof( uint16_type id,
                                uint16_type id_ent,
                                uint16_type node_id,
                                uint16_type entity_topodim,
                                uint16_type entity_id,
                                uint16_type entity_nsiblings,
                                uint16_type c = 0 )
        :
        M_id( id ),
        M_id_in_entity( id_ent ),
        M_node_id( node_id ),
        M_entity_topodim( entity_topodim ),
        M_entity_id( entity_id ),
        M_entity_nsiblings( entity_nsiblings ),
        M_c(c)
        {}

    FiniteElementDof& operator=( FiniteElementDof const& ) = default;
    FiniteElementDof& operator=( FiniteElementDof && ) = default;

    //! @return dof id
    constexpr uint16_type id() const noexcept { return M_id; }
    
    //! @return dof id in Entity
    constexpr uint16_type idInEntity() const noexcept { return M_id_in_entity; }
    
    //! @return node id associated to dof
    constexpr uint16_type nodeId() const noexcept { return M_node_id; }

    //! @return topological dimension of the entity associated to the dof
    constexpr uint16_type entityTopologicalDimension() const noexcept { return M_entity_topodim; }
    
    //! @return id of the entity associated to the dof
    constexpr uint16_type entityId() const noexcept { return M_entity_id; }

    //! @return number of dof siblings associated to the entity 
    constexpr uint16_type entityNumberOfSiblings() const noexcept { return M_entity_nsiblings; }
    
    //! @return dof component, default is 0
    constexpr uint16_type component() const noexcept { return M_c; }

private:
    uint16_type M_id, M_id_in_entity, M_node_id, M_entity_topodim, M_entity_id, M_entity_nsiblings, M_c;
    
};

inline std::ostream&
operator<<( std::ostream& os, FiniteElementDof const& d )
{
    os << "{" << d.id() << "," << d.idInEntity() << "," << d.entityTopologicalDimension() << "," << d.entityId() << "," << d.component() << "}";
    return os;
}
}
#endif
