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
#include <feel/feelcore/enumerate.hpp>
#include <feel/feelcore/serialization.hpp>
#include <feel/feeldiscr/functionspace.hpp>
#include <feel/feeldiscr/facet.hpp>

namespace Feel {

/**
 * @brief compute the set of operations @p ops of the element of a function space at DOF points over the range @p r
 * @ingroup Discretization
 *
 * @tparam Ts variadic template holding the type list of arguments
 *
 * @param _range range of element
 * @param _element element of a function space
 * @param _op std::vector of operations(lambda functions) to apply to the element
 *
 * @code {.cpp}
 * auto v = Xh->element();
 * auto res = minelt(_range=elements(mesh),_element=v,_op=[](auto const& a, auto const& b){ return a < b; });
 * // loop on results
 * for( auto item: res )
 * {
 *   auto [op_v, arg_op_v] = item;
 *   // do something
 * }
 * @endcode
 *
 * @return the vector of results and associated dof point coordinates
 */
template <typename... Ts>
auto opselt( Ts&&... v )
{
    auto args = NA::make_arguments( std::forward<Ts>( v )... );
    auto&& r = args.get( _range );
    auto&& e = args.get( _element );
    auto&& ops = args.get( _op ); // a vector of op

    using value_type = typename decay_type<decltype( e )>::value_type;
    using functionspace_type = typename decay_type<decltype( e )>::functionspace_type;
    using range_t = decay_type<decltype( r )>;
    using mesh_t = typename functionspace_type::mesh_type;
    constexpr int nRealDim = functionspace_type::nRealDim;
    value_type v_min = std::numeric_limits<value_type>::min();
    value_type v_max = std::numeric_limits<value_type>::max();
    std::vector<double> ops_v( ops.size() );
    for ( auto [i, op_v] : enumerate( ops_v ) )
        op_v = ops[i]( v_min, v_max ) ? v_max : v_min;
    std::vector<int> op_id( ops.size(), 0 );

    auto find_local_op = [&e, &ops, &ops_v, &op_id]( auto localdoftable, auto getDofId )
    {
        for ( auto const& ldof : localdoftable )
        {
            index_type dofpn = getDofId( ldof );
            value_type d = e( dofpn );

            for ( auto [i, op_v] : enumerate( ops_v ) )
            {
                if ( ops[i]( d, op_v ) )
                {
                    op_v = d;
                    op_id[i] = dofpn;
                }
            }
        }
    };

    for ( auto const& rangeElt : r )
    {
        auto const& meshElt = boost::unwrap_ref( rangeElt );

        if constexpr ( std::is_same_v<range_t, faces_reference_wrapper_t<mesh_t>> )
        {
            find_local_op( e.functionSpace()->dof()->faceLocalDof( meshElt.id() ), []( auto const& ldof )
                           { return ldof.index(); } );
        }
        else
        {
            find_local_op( e.functionSpace()->dof()->localDof( meshElt.id() ), []( auto const& ldof )
                           { return ldof.second.index(); } );
        }
    }
    std::vector<std::tuple<value_type, eigen_vector_type<nRealDim>>> res( ops.size() );
    for ( auto [i, op_v] : enumerate( ops_v ) )
    {
        int index = i;
        //auto [pt, thedof, comp] = e.functionSpace()->dof()->dofPoint( op_id[i] );
        auto pt = e.functionSpace()->dof()->dofPoint( op_id[i] ).template get<0>();
        eigen_vector_type<nRealDim> e_pt = emap<value_type>( pt );
        auto local_op = std::tuple{ op_v, e_pt };
        res[i] = mpi::all_reduce( e.functionSpace()->worldComm().globalComm(), local_op,
                                  [&ops, index]( auto const& a, auto const& b )
                                  { return ops[index]( std::get<0>( a ), std::get<0>( b ) ) ? a : b; } );
    }
    return res;
}


/**
 * @brief compute the min of the element of a function space at DOF points over the range @p r
 * @ingroup Discretization
 * 
 * @tparam Ts variadic template holding the type list of arguments
 *
 * @param _range range of element
 * @param _element element of a function space
 * @param _op operation to apply to the element
 *
 * @code {.cpp}
 * auto v = Xh->element();
 * auto [op_v,arg_op_v] = minelt(_range=elements(mesh),_element=v,_op=[](auto const& a, auto const& b){ return a < b; });
 * @endcode
 *
 * @return the value and the dof point coordinates
 */
template <typename... Ts>
auto opelt( Ts&&... v )
{
    auto args = NA::make_arguments( std::forward<Ts>( v )... );
    auto && r = args.get(_range);
    auto && e = args.get(_element);
    auto && op = args.get(_op); // no vector just a lambda function
    using value_type = typename decay_type<decltype( e )>::value_type;
    auto res = opselt( _range=r, _element=e, _op=std::vector<bool (*)(value_type const&,value_type const&)>{ { op } } );

    return res[0];
}
/**
 * @brief compute the min of the element of a function space at DOF points over the range @p r
 * @ingroup Discretization
 * 
 * @tparam RangeT range type to iterate on
 * @tparam ElementT type of finite element function
 *
 * @param r range of element
 * @param e element of a function space
 *
 * @code {.cpp}
 * auto v = Xh->element();
 * auto [min_v,arg_min_v] = minelt(_range=elements(mesh),_element=v);
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
    return  opelt( _range=r, _element=e, _op=[]( auto const& a, auto const& b ){ return a < b; } );
}
/**
 * @brief compute the max of the element of a function space at DOF points over the range @p r
 * @ingroup Discretization
 * 
 * @tparam RangeT range type to iterate on
 * @tparam ElementT type of finite element function
 *
 * @param r range of element
 * @param e element of a function space
 *
 * @code {.cpp}
 * auto v = Xh->element();
 * auto [max_v,arg_max_v] = maxelt(_range=elements(mesh),_element=v);
 * @endcode
 *
 * @return the max value and the dof point coordinates
 */
template <typename... Ts>
auto maxelt( Ts&&... v )
{
    auto args = NA::make_arguments( std::forward<Ts>( v )... );
    auto && r = args.get(_range);
    auto && e = args.get(_element);
    return  opelt( _range=r, _element=e, _op=[]( auto const& a, auto const& b ){ return a > b; } );
}


/**
 * @brief compute the min of the element of a function space at DOF points over the range @p r
 * @ingroup Discretization
 *
 * @tparam RangeT range type to iterate on
 * @tparam ElementT type of finite element function
 *
 * @param r range of element
 * @param e element of a function space
 *
 * @code {.cpp}
 * auto v = Xh->element();
 * auto r = minmaxelt(_range=elements(mesh),_element=v);
 * for( auto const& [o_v,arg_o_v] : enumerate(r) )
 * {
 *      std::cout << "i: " << o_v << ", " << arg_o_v << std::endl;
 * }
 * r[0]; // min
 * r[1]; // max
 * @endcode
 *
 * @return the min value and the dof point coordinates
 */
template <typename... Ts>
auto minmaxelt( Ts&&... v )
{
    auto args = NA::make_arguments( std::forward<Ts>( v )... );
    auto && r = args.get(_range);
    auto && e = args.get(_element);
    using value_type = typename decay_type<decltype( e )>::value_type;
    return  opselt( _range=r, _element=e, _op=std::vector<bool (*)(value_type const&,value_type const&)>{ 
                                                            { []( value_type const& a, value_type const& b ){ return a < b; },
                                                              []( value_type const& a, value_type const& b ){ return a > b; } } } );
}
    
} // namespace Feel