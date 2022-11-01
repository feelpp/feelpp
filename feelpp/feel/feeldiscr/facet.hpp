/**
 * @file facet.hpp
 * @author Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 * @brief 
 * @version 0.1
 * @date 2022-07-21
 * 
 * @copyright Copyright (c) 2022 Université de Strasbourg
 * 
 */
#pragma once

namespace Feel {

/**
 * @brief get the local indexing of a facet from the global id
 * @ingroup Discretization
 * @tparam FacetType type of face
 * @param f facet 
 * @return a tuple containing the element 0 or 1 to which the facet belongs and the local index of the facet in this element
 */
template<typename FacetType>
auto facetGlobalToLocal( FacetType && f )
{
    int face_id = f.pos_first();
    int faceConnectionId = 0;
    
    if ( f.element( 0 ).isGhostCell() )
    {
        face_id = std::forward<FacetType>(f).pos_second();
        faceConnectionId = 1;
    }
    return std::tuple{faceConnectionId, face_id};
}

/**
 * @brief get the local indexing of a facet from the global id and a doftable
 * @ingroup Discretization
 * 
 * @warning if the tuple {-1,-1} is returned it means with are on a partial mesh and the facet is in it
 * 
 * @tparam FacetType type of facet
 * @tparam DofTableType type of doftable
 * @param f facet
 * @param dof doftable
 * @return a tuple containing the element 0 or 1 to which the facet belongs and the local index of the facet in this element
 */
template <typename FacetType, typename DofTableType>
auto facetGlobalToLocal( FacetType&& f, DofTableType&& dof )
{
    bool hasMeshSupportPartial = dof->hasMeshSupport() && dof->meshSupport()->isPartialSupport();
    bool hasDofTableMPIExtended = dof->buildDofTableMPIExtended();

    int face_id = f.pos_first();
    int faceConnectionId = 0;

    if ( hasMeshSupportPartial )
    {
        auto const& elt0 = f.element( 0 );
        if ( !dof->meshSupport()->hasElement( elt0.id() ) || ( !hasDofTableMPIExtended && elt0.isGhostCell() ) )
        {
            if ( !f.isConnectedTo1() )
                return std::tuple{-1, -1};
            auto const& elt1 = f.element( 1 );
            if ( !dof->meshSupport()->hasElement( elt1.id() ) || ( !hasDofTableMPIExtended && elt1.isGhostCell() ) )
                return std::tuple{-1, -1};
            face_id = f.pos_second();
            faceConnectionId = 1;
        }
    }
    else if ( !hasDofTableMPIExtended && f.element( 0 ).isGhostCell() )
    {
        face_id = std::forward<FacetType>( f ).pos_second();
        faceConnectionId = 1;
    }
    return std::tuple{ faceConnectionId, face_id };
}
} // Feel