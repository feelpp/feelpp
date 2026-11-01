/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileContributor: Vincent Chabannes <vincent.chabannes@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2014-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
   \file integrate.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2014-01-04
 */
#ifndef FEELPP_VF_INTEGRATE_H
#define FEELPP_VF_INTEGRATE_H

#include <feel/feelvf/expr.hpp>
#include <feel/feelvf/integrator.hpp>

namespace Feel {

namespace detail
{
template<typename ArgExprType,typename ArgRangeType,typename ArgQuadType,typename ArgQuad1Type>
struct integrate_type
 {
    using _expr_type = std::decay_t<ArgExprType>;
    using _range_type = typename Feel::detail::quadptlocrangetype< std::decay_t<ArgRangeType> >::type;

    using _element_iterator = typename _range_type::iterator_t;
    using _element_type = typename _range_type::iterator_t;
    static inline const uint16_type geoOrder = _element_type::nOrder;
    
    using expr_order_t = Feel::vf::ExpressionOrder<_range_type,_expr_type>;
    //using _value_type = typename _expr_type::value_type;
    using im_default_type = im_t<typename expr_order_t::the_element_type, typename _expr_type::value_type>;
    using __quad_type = std::decay_t<ArgQuadType>;
    using __quad1_type = std::decay_t<ArgQuad1Type>;
    using _im_type = Feel::vf::detail::integrate_im_type<_range_type,_expr_type,__quad_type,__quad1_type>;
    using _quad_type = typename _im_type::_quad_type;
    using _quad1_type = typename _im_type::_quad1_type;
    typedef Expr<Integrator<_range_type, _quad_type, _expr_type, _quad1_type> > expr_type;
    typedef std::shared_ptr<QuadPtLocalization<_range_type,_quad_type,_expr_type > > _quadptloc_ptrtype;
 };

} // detail


template <typename ... Ts>
auto integrate( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    auto && range = args.get(_range);
    auto && expr = args.get(_expr);
    auto && quad = args.get_else(_quad,quad_order_from_expression );
    auto && quad1 = args.get_else(_quad1,quad_order_from_expression );
    GeomapStrategyType geomap = args.get_else(_geomap,GeomapStrategyType::GEOMAP_OPT);
    auto geometryEvaluation = args.get_else(
        _geometry_evaluation,
        geometry::geometry_evaluation_strategy::automatic );
    bool use_harts = args.get_else(_use_harts,false);
    auto localFormLowering = args.get_else(
        _local_form_lowering, vf::local_form_lowering::automatic );
    bool verbose = args.get_else(_verbose,false);

    using _integrate_helper_type = Feel::detail::integrate_type<decltype(expr),decltype(range),decltype(quad),decltype(quad1)>;
    auto && quadptloc = args.get_else(_quadptloc, typename _integrate_helper_type::_quadptloc_ptrtype{});

    auto the_ims = _integrate_helper_type::_im_type::im( quad,quad1,expr );
    auto const& the_im = the_ims.first;
    auto const& the_im1 = the_ims.second;

    auto ret = integrate_impl( range, the_im, expr,
                               Feel::detail::geomapStrategy( range, geomap ),
                               geometryEvaluation, the_im1,
                               use_harts, quadptloc, localFormLowering );

    if ( verbose )
    {
        std::cout << " -- integrate: size(range) = " << std::distance( ret.expression().beginElement(),
                                                                       ret.expression().endElement() ) << "\n";
        std::cout << " -- integrate: quad = " << ret.expression().im().nPoints() << "\n";
        std::cout << " -- integrate: quad1 = " << ret.expression().im2().nPoints() << "\n";
        //std::cout << " -- integrate: geomap = " << geomap << "\n";
    }

    return ret;
}

} // namespace Feel

#endif
