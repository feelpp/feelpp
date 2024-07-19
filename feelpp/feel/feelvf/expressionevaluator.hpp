/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-11-24

  Copyright (C) 2009-2012 Universite Joseph Fourier (Grenoble I)
  Copyright (C) 2011-present Feel++ Consortium
  
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

#ifndef FEELPP_EXPRESSION_EVALUATOR_HPP
#define FEELPP_EXPRESSION_EVALUATOR_HPP

#include <feel/feelpoly/context.hpp>
#include <feel/feelpoly/geomap.hpp>

namespace Feel
{
template<typename RangeType>
class ExpressionEvaluatorBase
{
public:
    using range_type = RangeType;
    using element_type = typename range_type::element_t;
    using gm_type = typename element_type::gm_type;
    using parameterspace_type = ParameterSpace<>;
    using parameterspace_ptrtype = std::shared_ptr<parameterspace_type>;
    using parameterelement_type = typename parameterspace_type::element_type;

public:
    ExpressionEvaluatorBase(range_type const& r, int comp = 0) : M_range(r), M_comp(comp) {};
    virtual void init(int o) = 0;
    virtual void update(element_type const& elt) = 0;
    virtual bool update(parameterelement_type const& mu) = 0;
    virtual int nPoints() = 0;
    virtual double weight(int i) = 0;
    virtual double eval(int q, int comp) = 0;
    virtual int order() = 0;
    int component() { return M_comp; }

protected:
    range_type M_range;
    int M_comp;
};

template<typename EltT, typename ExprT>
class ExpressionEvaluator : public ExpressionEvaluatorBase<EltT>
{
protected:
    using super = ExpressionEvaluatorBase<EltT>;
    using range_type = typename super::range_type;
    using element_type = typename super::element_type;
    using gm_type = typename super::gm_type;

    using expr_type = ExprT;
    using context_type = typename gm_type::template Context<std::remove_const_t<element_type>>;
    using context_ptrtype = std::shared_ptr<context_type>;
    using map_gmc_type = Feel::vf::map_gmc_type<context_type>;
    using evaluator_type = typename expr_type::template tensor<map_gmc_type>;
    using evaluator_ptrtype = std::shared_ptr<evaluator_type>;
    using weights_type = ublas::vector<double>;

    using parameterspace_type = typename super::parameterspace_type;
    using parameterspace_ptrtype = typename super::parameterspace_ptrtype;
    using parameterelement_type = typename super::parameterelement_type;

public:
    ExpressionEvaluator( range_type const& r, expr_type const& ex, int comp = 0);
    void init(int o) override;
    void update(element_type const& elt) override;
    virtual bool update(parameterelement_type const& mu) override { return true; }
    int nPoints() override;
    double weight(int i) override;
    double eval(int q, int comp) override;
    int order() override;
    ExprT expression() const { return M_expr; }

    expr_type M_expr;

protected:
    context_ptrtype M_ctx;
    evaluator_ptrtype M_evaluator;
    weights_type M_weights;
};

template<typename RangeType, typename ExprT>
ExpressionEvaluator<RangeType, ExprT>::ExpressionEvaluator( range_type const& r, expr_type const& ex, int comp )
    : super(r, comp),
      M_expr(ex)
{}

template<typename RangeType, typename ExprT>
void
ExpressionEvaluator<RangeType, ExprT>::init(int o)
{
    auto const& eltForInit = this->M_range.front();
    auto q = _Q(o);

    using iim_type = vf::detail::integrate_im_type<decltype(this->M_range),decltype(M_expr),decltype(q),decltype(q)>;
    auto ims = iim_type::im(q,q,M_expr);
    auto im = ims.first;
    auto pts = im.points();
    auto gm = eltForInit.gm();
    auto geopc = gm->preCompute( pts );

    M_ctx = gm->template context<vm::POINT|vm::JACOBIAN|expr_type::context>( eltForInit, geopc );
    M_evaluator = std::make_shared<evaluator_type>(M_expr.evaluator( mapgmc(M_ctx) ) );
    M_weights = im.weights();
}

template<typename RangeType, typename ExprT>
void
ExpressionEvaluator<RangeType, ExprT>::update(element_type const& eltWrap)
{
    auto const& elt = unwrap_ref( eltWrap );
    M_ctx->template update<vm::POINT | vm::JACOBIAN | expr_type::context>( elt );
    M_evaluator->update( vf::mapgmc( M_ctx ) );
}

template<typename RangeType, typename ExprT>
double
ExpressionEvaluator<RangeType, ExprT>::eval(int q, int comp)
{
    return M_evaluator->evalq(comp,0,q)*M_ctx->J(q);
}

template<typename RangeType, typename ExprT>
int
ExpressionEvaluator<RangeType, ExprT>::nPoints()
{
    return M_ctx->nPoints();
}

template<typename RangeType, typename ExprT>
double
ExpressionEvaluator<RangeType, ExprT>::weight(int i)
{
    return M_weights[i];
}

template<typename RangeType, typename ExprT>
int
ExpressionEvaluator<RangeType, ExprT>::order()
{
    return M_expr.polynomialOrder();
}

template<typename EltT, typename ExprT>
class ExpressionEvaluatorParam : public ExpressionEvaluator<EltT, ExprT>
{
protected:
    using super = ExpressionEvaluator<EltT,ExprT>;
    using range_type = typename super::range_type;
    using element_type = typename super::element_type;
    using gm_type = typename super::gm_type;

    using expr_type = ExprT;
    using context_type = typename super::context_type;
    using context_ptrtype = typename super::context_ptrtype;
    using map_gmc_type = typename super::map_gmc_type;
    using evaluator_type = typename super::evaluator_type;
    using evaluator_ptrtype = typename super::evaluator_ptrtype;
    using weights_type = typename super::weights_type;

    using parameterspace_type = typename super::parameterspace_type;
    using parameterspace_ptrtype = typename super::parameterspace_ptrtype;
    using parameterelement_type = typename super::parameterelement_type;

    using ExpressionEvaluator<EltT, ExprT>::update;

public:
    ExpressionEvaluatorParam( range_type const& r, expr_type const& ex,
                              parameterelement_type const& mu, int comp = 0 ):
        super(r,ex,comp),
        M_mu(mu)
        {}

    virtual bool update(parameterelement_type const& mu) override;

protected:
    parameterelement_type M_mu;
};

template<typename EltT, typename ExprT>
bool ExpressionEvaluatorParam<EltT, ExprT>::update(parameterelement_type const& mu)
{
    M_mu = mu;
    return true;
}

template<typename EltT, typename ExprT, typename FctT>
class ExpressionEvaluatorNonLinear : public ExpressionEvaluatorParam<EltT, ExprT>
{
public:
    using super = ExpressionEvaluatorParam<EltT,ExprT>;
    using range_type = typename super::range_type;
    using element_type = typename super::element_type;
    using gm_type = typename super::gm_type;

    using expr_type = ExprT;
    using context_type = typename super::context_type;
    using context_ptrtype = typename super::context_ptrtype;
    using map_gmc_type = typename super::map_gmc_type;
    using evaluator_type = typename super::evaluator_type;
    using evaluator_ptrtype = typename super::evaluator_ptrtype;
    using weights_type = typename super::weights_type;

    using parameterspace_type = typename super::parameterspace_type;
    using parameterspace_ptrtype = typename super::parameterspace_ptrtype;
    using parameterelement_type = typename super::parameterelement_type;

    using function_element_type = FctT;
    using update_function_type = std::function<bool(parameterelement_type const&, function_element_type&)>;

    using ExpressionEvaluator<EltT, ExprT>::update;

public:
    ExpressionEvaluatorNonLinear( range_type const& r, expr_type const& ex,
                                  parameterelement_type& mu, FctT& u, int comp = 0 )
        : super(r,ex,mu,comp),
              M_u(u)
        {}

    bool update(parameterelement_type const& mu) override;

    update_function_type M_fct;
protected:
    function_element_type& M_u;
};

template<typename EltT, typename ExprT, typename FctT>
bool
ExpressionEvaluatorNonLinear< EltT, ExprT, FctT>::update( parameterelement_type const& mu )
{
    this->M_mu = mu;
    if( M_fct )
    {
        if( M_fct(this->M_mu, M_u) )
            this->M_expr(idv(M_u));
        else
            return false;
    }
    return true;
}
} // Feel

#endif
