/**
 * @file context.hpp
 * @author Christophe Prud'homme (christophe.prudhomme@cemosis.fr)
 * @brief 
 * @version 0.1
 * @date 2022-07-24
 * 
 * @copyright Copyright (c) 2022  Feel++ Consortium
 * @copyright Copyright (c) 2022 Université de Strasbourg
 * 
 */
#pragma once


#include <feel/feelpoly/context.hpp>
#include <feel/feelmesh/enums.hpp>

namespace Feel {

/**
 * @brief precompute geometric mapping at a set of points
 * 
 */
template <typename... Ts>
[[nodiscard]] auto geomapPrecompute( Ts&&... v ) 
{
    auto args = NA::make_arguments( std::forward<Ts>( v )... );
    auto&& element = args.get( _element );
    auto&& type = args.get( _type );
    auto&& geomap = args.get( _geomap );
    auto&& pointset = args.get( _pointset ); // lambda function

    using element_t = decay_type<decltype( element )>;
    using gm_t = decay_type<decltype( geomap )>;
    using geopc_t = typename gm_t::precompute_type;
    using geopc_ptr_t = std::shared_ptr<geopc_t>;
    
    if constexpr(  std::is_same_v<decay_type<decltype(type)>, on_facets_t> )
    {
        using permutation_t = typename element_t::permutation_type;//typename element_t::template PermutationSubEntity<1>;

        std::vector<std::map<permutation_t, geopc_ptr_t>> geopc( element_t::numTopologicalFaces );

        for ( uint16_type __f = 0; __f < element_t::numTopologicalFaces; ++__f )
        {
            for ( permutation_t __p( permutation_t::IDENTITY );
                    __p < permutation_t( permutation_t::N_PERMUTATIONS ); ++__p )
            {
                geopc[__f][__p] = std::make_shared<geopc_t>( geomap, pointset( __f, __p ) );
            }
        }
        return geopc;
    }
    else
    {
        return std::make_shared<geopc_t>( geomap, pointset(0,0) );
    }
}

/**
 * @brief get the the geometric mapping context associated to a functionspace
 * 
 * @tparam Ts argument type variadic list 
 * @param v argument variadic list
 * 
 * @code
 * auto Xh=Pch<2>(mesh);
 * auto ctx = context(_space=Xh,_type=on_elements_t());
 * @endcode
 * 
 * @return the context
 */
template <typename... Ts>
[[nodiscard]] auto context( Ts&&... v )
{
    auto args = NA::make_arguments( std::forward<Ts>( v )... );
    auto&& space = args.get( _space );
    auto&& type = args.get( _type );

    auto&& pointset = args.get_else( _pointset, [&space]( int, int) { return space->fe()->points(); } );
    auto&& context_at_compile_time = args.get_else( _context_at_compile_time, POINT );
    auto&& context_at_run_time = args.get_else( _context, 0 );

    auto mesh = space->mesh();
    using gm_type = typename decay_type<decltype( mesh )>::gm_type;
    if constexpr ( std::is_same_v<decay_type<decltype(type)>, on_elements_t> ) 
    {
        auto const& initElt = mesh->beginElement()->second;
        auto geopc = geomapPrecompute( _element=initElt, _type=on_elements_t(), _geomap=space->gm(), _pointset=pointset );
        auto ctx = space->gm()->template context<context_at_compile_time>( space->gm(), initElt, geopc, context_at_run_time );
        return ctx;
    }
    else if constexpr ( std::is_same_v<decay_type<decltype( type )>, on_facets_t> )
    {
        using geoelement_t = typename decay_type<decltype(mesh)>::element_type;
        auto const& initElt = mesh->beginElement()->second;
        auto const& initFace = mesh->beginFace()->second;
        auto geopc = geomapPrecompute(
            _element = initElt, _type = on_facets_t(), _geomap = space->gm(), _pointset = pointset );
  
        int face_id = initFace.pos_first();
     
        auto ctx = space->gm()->template context<context_at_compile_time>( space->gm(), initFace.element( 0 ), geopc, face_id, context_at_run_time );

        return ctx;        
    }
}

} // Feel