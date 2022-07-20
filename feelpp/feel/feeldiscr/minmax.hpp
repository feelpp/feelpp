/**
 * @file minmax.hpp
 * @author Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 * @brief min and max functions of an element of a function space
 * @version 0.1
 * @date 2022-06-21
 * 
 * @copyright Copyright (c) 2022 Feel++ Consortium
 * @copyright Copyright (c) 2022 Université de Strasbourg
 * 
 */
#pragma once

#include <feel/feelalg/glas.hpp>
#include <feel/feelcore/parameter.hpp>
#include <feel/feelcore/serialization.hpp>
#include <feel/feeldiscr/functionspace.hpp>

namespace Feel {

/**
 * @brief compute the min of the element of a function space at DOF points over the range @p r
 * 
 * @tparam RangeT range type to iterate on
 * @tparam ElementT type of finite element function
 * 
 * @param r range of element
 * @param e element of a function space
 * 
 * @code {.cpp}
 * auto v = Xh->element();
 * auto [min_v,arg_min_v] = min(_range=elements(mesh),_element=v);
 * @endcode
 * 
 * @return the min value and the dof point coordinates
 */
template <typename... Ts>
auto minelt( Ts&&... v )
{
    auto args = NA::make_arguments( std::forward<Ts>( v )... );
    auto && r = args.get(_range);
    auto && e = args.get(_element);

    using value_type = typename decay_type<decltype(e)>::value_type;
    using functionspace_type = typename decay_type<decltype(e)>::functionspace_type;
    constexpr int nRealDim = functionspace_type::nRealDim;
    double min_v{std::numeric_limits<value_type>::max()};
    int min_id{ 0 };

    for ( auto const& rangeElt : r )
    {
        auto const& meshElt = boost::unwrap_ref( rangeElt );
        index_type eid = meshElt.id();

        for ( uint16_type local_id=0; local_id < functionspace_type::fe_type::nLocalDof; ++local_id )
        {
            int dofpn = e.functionSpace()->dof()->localToGlobal( eid, local_id, 0 ).index();
            value_type d =e( dofpn );
            if ( min_v > d )
            {
                min_v = d;
                min_id = dofpn;
            }
        }
    }
    auto [pt,thedof,comp] = e.functionSpace()->dof()->dofPoint(min_id);
    eigen_vector_type<nRealDim> e_pt = emap<value_type>( pt );
    auto local_min = std::tuple{min_v,e_pt};
    std::vector<std::tuple<value_type,eigen_vector_type<nRealDim>>> all_min(e.functionSpace()->worldComm().globalSize());
    mpi::all_gather( e.functionSpace()->worldComm().globalComm(), local_min, all_min );
    return *std::min_element( all_min.begin(), all_min.end(), 
                                []( auto const& a, auto const& b  ){
                                    return std::get<0>(a) < std::get<0>(b); }
                             );
}

/**
 * @brief compute the min of the element over the range @p r
 * 
 * @tparam RangeT range type to iterate on
 * @return the min value and the dof point coordinates
 */
template <typename RangeT, typename ElementT>
std::tuple<typename ElementT::value_type, typename ElementT::value_type, eigen_vector_type<ElementT::nRealDim>, eigen_vector_type<ElementT::nRealDim>> minmax( RangeT const& r )
{
#if 0
    using value_type = typename ElementT::value_type;
    eigen_matrix_type<2> minmax_v{std::numeric_limits<value_type>::max(),std::numeric_limits<value_type>::min()};
    std::pair<int,int> minmax_id{ 0, 0 };
    auto& [min_v_id,max_v_id] = minmax_id;
    for ( auto const& rangeElt : r )
    {
        auto const& meshElt = boost::unwrap_ref( rangeElt );
        index_type eid = meshElt.id();

        for ( uint16_type local_id=0; local_id < functionspace_type::fe_type::nLocalDof; ++local_id )
        {
            int dofpn = e.functionSpace()->dof()->localToGlobal( eid, local_id, 0 ).index();
            value_type d = this->operator()( dofpn );
            if ( minmax_v[0] > d )
            {
                minmax_v[0] = d;
                min_v_id = dofpn;
            }
            if ( minmax_v[1] < d )
            {
                minmax_v[1] = d;
                max_v_id = dofpn;
            }
        }
    }
    auto local_minmax = std::tuple{minmax_v,e.functionspace()->dof()->dofPoint(min_id).template get<0>()};
    std::vector<std::tuple<value_type,eigen_vector_type<nRealDim>>> all_minmax;
    mpi::all_gather( e.functionSpace()->worldComm(),
                        local_minmax,
                        all_minmax );
    return std::pair{ *std::min_element( all_minmax.begin(), all_minmax.end(), 
                                []( auto const& a, auto const& b  ){
                                    return a.get<0>()[0] < b.get<0>()[0];
                            } ),
                        *std::max_element( all_minmax.begin(), all_minmax.end(), 
                                []( auto const& a, auto const& b  ){
                                    return a.get<0>()[1] < b.get<0>()[1];
                            } ) };
#endif                            
}
}