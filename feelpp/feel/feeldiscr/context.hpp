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
 * auto ctx = context(_geomap=Xh->gm(),_type=on_elements_t(),_pointset=Xh->fe()->points());
 * @endcode
 * 
 * @return the context
 */
template <typename... Ts>
[[nodiscard]] auto context( Ts&&... v )
{
    auto args = NA::make_arguments( std::forward<Ts>( v )... );
    auto&& geomap = args.get( _geomap );
    auto&& type = args.get( _type );
    auto&& element = args.get( _element );
    auto&& pointset = args.get( _pointset );
    constexpr size_type context_at_compile_time = args.get_else( _context_at_compile_time, POINT|NORMAL|JACOBIAN );
    auto&& context_at_run_time = args.get_else( _context, 0 );
    using gm_type = decay_type<decltype( geomap )>;

    
    if constexpr ( std::is_same_v<decay_type<decltype(type)>, on_elements_t> ) 
    {
        auto geopc = geomapPrecompute( _element=element, _type=on_elements_t(), _geomap=geomap, _pointset=pointset );
        auto ctx = geomap->template context<context_at_compile_time>( geomap, element, geopc, context_at_run_time );
        return ctx;
    }
    else if constexpr ( std::is_same_v<decay_type<decltype( type )>, on_facets_t> )
    {
        auto geopc = geomapPrecompute(
            _element = element.element(0), _type = on_facets_t(), _geomap = geomap, _pointset = pointset );
  
        int face_id = element.pos_first();
     
        auto ctx = geomap->template context<context_at_compile_time>( geomap, element.element( 0 ), geopc, face_id, context_at_run_time );

        return ctx;        
    }
}

} // Feel