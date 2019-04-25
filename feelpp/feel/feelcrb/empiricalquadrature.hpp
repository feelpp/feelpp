/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-11-24

  Copyright (C) 2009-2012 Universite Joseph Fourier (Grenoble I)

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

#ifndef FEELPP_EQ_HPP
#define FEELPP_EQ_HPP

#include <feel/feelcrb/parameterspace.hpp>

#if defined(FEELPP_HAS_GLPK_H)
#include <glpk.h>
#endif /* FEELPP_HAS_GLPK_H */

using namespace Feel;

template<typename RangeType>
class ExpressionEvaluatorBase
{
protected:
    using range_type = RangeType;
    using element_type = typename boost::tuples::element<1_c,range_type>::type::value_type;
    using gm_type = typename element_type::type::gm_type;
public:
    ExpressionEvaluatorBase(range_type const& r) : M_range(r) {};
    virtual void init(int o) = 0;
    virtual void update(element_type elt) = 0;
    virtual int nPoints() = 0;
    virtual double weight(int i) = 0;
    virtual double eval(int q, int comp) = 0;
    virtual int order() = 0;

protected:
    range_type M_range;
};

template<typename EltT, typename ExprT>
class ExpressionEvaluator : public ExpressionEvaluatorBase<EltT>
{
    using super = ExpressionEvaluatorBase<EltT>;
    using range_type = typename super::range_type;
    using element_type = typename super::element_type;
    using gm_type = typename super::gm_type;

    using expr_type = ExprT;
    using context_type = typename gm_type::template Context<vm::POINT|vm::JACOBIAN|expr_type::context,typename std::remove_const<typename element_type::type>::type >;
    using context_ptrtype = std::shared_ptr<context_type>;
    using map_gmc_type = map_gmc_type<context_type>;
    using evaluator_type = typename expr_type::template tensor<map_gmc_type>;
    using evaluator_ptrtype = std::shared_ptr<evaluator_type>;
    using weights_type = ublas::vector<double>;

public:
    ExpressionEvaluator( range_type const& r, expr_type const& ex );
    void init(int o) override;
    void update(element_type elt) override;
    int nPoints() override;
    double weight(int i) override;
    double eval(int q, int comp) override;
    int order() override;

    expr_type M_expr;

private:
    context_ptrtype M_ctx;
    evaluator_ptrtype M_evaluator;
    weights_type M_weights;
};

template<typename RangeType>
class EmpiricalQuadrature
{
    using range_type = RangeType;
    using element_type = typename boost::tuples::element<1_c,range_type>::type::value_type;
    using gm_type = typename element_type::type::gm_type;

    using expressionevalbase_type = ExpressionEvaluatorBase<element_type>;
    using expressionevalbase_ptrtype = std::shared_ptr<expressionevalbase_type>;
    using expressionevalbase_vectype = std::vector<expressionevalbase_ptrtype>;
    template<typename ExprT>
    using expressioneval_type = ExpressionEvaluator<element_type, ExprT>;
    // using exprs_type = std::tuple<ExprTs...>;
    // using contexts_type = std::tuple<typename gm_type::template Context<vm::POINT|vm::JACOBIAN|ExprTs::context,typename std::remove_const<typename element_type::type>::type >...>;
    // using contexts_ptrtype = std::tuple<std::shared_ptr<typename gm_type::template Context<vm::POINT|vm::JACOBIAN|ExprTs::context,typename std::remove_const<typename element_type::type>::type >>...>;
    // // using map_gmcs_type = std::tuple<map_gmc_type<typename gm_type::template Context<vm::POINT|vm::JACOBIAN|ExprTs::context,typename std::remove_const<typename element_type::type>::type >>...>;
    // // using evaluators_type = typename ExprTs::template tensor<map_gmc_type>...;
    // // using evaluators_ptrtype = std::shared_ptr<evaluators_type...>;
    // // using evaluators_tuple_type = std::tuple<evaluators_ptrtype...>;

    // using expr_type = typename std::tuple_element<0,exprs_type>::type;
    // using context_type = typename gm_type::template Context<vm::POINT|vm::JACOBIAN|expr_type::context,typename std::remove_const<typename element_type::type>::type >;
    // using context_ptrtype = std::shared_ptr<context_type>;
    // using map_gmc_type = map_gmc_type<context_type>;
    // using evaluator_type = typename expr_type::template tensor<map_gmc_type>;
    // using evaluator_ptrtype = std::shared_ptr<evaluator_type>;

    using parameterspace_type = ParameterSpace<>;
    using parameterspace_ptrtype = std::shared_ptr<parameterspace_type>;
    using parameterelement_type = parameterspace_type::element_type;
    using sampling_type = typename parameterspace_type::sampling_type;
    using sampling_ptrtype = typename parameterspace_type::sampling_ptrtype;

    using evaluation_type = std::vector<std::vector<std::vector<double> > >;
    using local_result_type = std::vector<std::vector<double> >;

public:
    EmpiricalQuadrature( range_type const& range, parameterelement_type& mu,
                         parameterspace_ptrtype const& Dmu,
                         int comp = 0, std::string const& prefix = "" );
    template<typename ExprT>
    void addExpression(ExprT const& ex);
    int offline();
    double evaluate( int m = 0);

    int dimension() const { return M_weights.size(); }

private:
    std::vector<int> solveLP( evaluation_type const& eval, local_result_type const& res );

private:
    range_type M_range;
    expressioneval_vectype M_exprevals;
    // exprs_type M_exprs;
    // contexts_ptrtype M_ctxs;
    // expr_type M_expr;
    // context_ptrtype M_ctx;
    // evaluator_ptrtype M_evaluator;
    int M_component;
    std::vector<std::vector<double> > M_evaluations;
    std::vector<double> M_weights;
    std::vector<std::pair<int,element_type> >  M_points;
    parameterspace_ptrtype M_Dmu;
    sampling_ptrtype M_trainset;
    parameterelement_type M_mu;

    std::string M_prefix;
    double M_tol;
    int M_M;
    int M_J;
    int M_numElts;
    int M_N;
};

template<typename RangeType>
EmpiricalQuadrature<RangeType>::EmpiricalQuadrature( range_type const& range,
                                                     parameterelement_type& mu,
                                                     parameterspace_ptrtype const& Dmu,
                                                     int comp, std::string const& prefix )
    :
    M_range(range),
    // M_exprs(exprs...),
    // M_expr(std::get<0>(M_exprs)),
    M_component(comp),
    M_Dmu(Dmu),
    M_mu(mu),
    M_prefix(prefix),
    M_tol(doption(_prefix=M_prefix,_name="eq.tolerance")),
    M_M(0),
    M_J(ioption(_prefix=M_prefix,_name="eq.sampling-size"))
{
    M_numElts = nelements(M_range);
    M_trainset = M_Dmu->sampling();
    M_trainset->randomize(M_J);
}

template<typename RangeType>
template<typename ExprT>
void
EmpiricalQuadrature<RangeType>::addExpression( ExprT const& ex)
{
    auto exev = std::make_shared<expressioneval_type>(ex);
    M_exprevals.push_back(exev);
    M_M++;
}

template<typename RangeType>
int
EmpiricalQuadrature<RangeType>::offline()
{
    if( M_numElts <= 0 )
        return 1;

    auto const eltForInit = boost::unwrap_ref(*boost::get<1>(M_range));

    int max_order = std::reduce(M_exprevals.first, M_exprevals.end(), 0.0,
                                [](auto const& a, auto const& b)
                                    {
                                        return std::max(a.order(),b.order());
                                    });
    for( auto const& ee : M_exprevals )
        ee.init(max_order);
    // std::apply([max_order](auto const& ...ex) mutable {
    //                return (max_order = ... = std::max(max_order, (int)ex.polynomialOrder()));},
    //            M_exprs);
    // auto q = _Q(max_order);

    // using iim_type = vf::detail::integrate_im_type<decltype(M_range),decltype(M_expr),decltype(q),decltype(q)>;
    // auto ims = iim_type::im(q,q,M_expr);
    // auto im = ims.first;
    // auto pts = im.points();
    // int npts = im.nPoints();
    // auto gm = eltForInit.gm();
    // auto geopc = gm->preCompute( pts );

    // M_ctxs = std::make_tuple(gm->context( eltForInit, geopc ));
    // M_ctx = gm->template context<vm::POINT|vm::JACOBIAN|expr_type::context>( eltForInit, geopc );
    // M_evaluator = std::make_shared<evaluator_type>(M_expr.evaluator( mapgmc(M_ctx) ) );

    int nPts = M_exprevals.begin()->nPoints();
    M_N = M_numElts*nPts;
    evaluation_type eval(M_M, std::vector<std::vector<double> >(M_J, std::vector<double>(M_N, 0.0)));
    local_result_type res(M_M, std::vector<double>(M_N, 0.0) );
    int m = 0, j = 0, n = 0;
    for( auto const& eltWrap : M_range )
    {
        auto const& elt = unwrap_ref( eltWrap );
        if ( elt.processId() != Environment::rank() )
            continue;

        for( auto const& ee : M_exprevals )
            ee.update( elt );
        // M_ctx->update( elt );
        // M_evaluator->update( vf::mapgmc( M_ctx ) );

        for ( uint16_type q = 0; q < nPts; ++q )
        {
            for( m = 0; m < M_M; ++m )
            {
                for( j = 0; j < M_J; ++j )
                {
                    M_mu = M_trainset->at(j);
                    eval[m][j][n] = M_exprevals[m]->eval(q, M_component); //M_evaluator->evalq( M_component,0,q )*M_ctx->J(q);
                    res[m][j] += M_exprevals[m]->weight(q)*eval[m][j][n];
                }
            }
            n++;
        }
    }

    // double b = 0.0;
    // mpi::all_reduce(Environment::worldComm().globalComm(), res[0][0], b, std::plus<double>());
    // Feel::cout << "res = " << b << std::endl;

    auto indexes = solveLP(eval,res);

    Feel::cout << "non zero : " << indexes.size() << std::endl;
    if( indexes.size() == 0 )
        return 2;

    n = 0;
    auto it = indexes.begin();
    for( auto const& eltWrap : M_range )
    {
        auto const& elt = unwrap_ref( eltWrap );
        if ( elt.processId() != Environment::rank() )
            continue;
        // M_ctx->update( elt );
        for ( uint16_type q = 0; q < nPts; ++q )
        {
            if( it != indexes.end() && n == *it )
            {
                ++it;
                M_points.push_back(std::make_pair(q,eltWrap));
            }
            n++;
        }
    }

    return 0;
}

template<typename RangeType>
std::vector<int>
EmpiricalQuadrature<RangeType>::solveLP( evaluation_type const& eval, local_result_type const& res )
{
    tic();
    glp_prob *lp;
    int ia[1+M_M*M_J*M_N], ja[1+M_M*M_J*M_N];
    double ar[1+M_M*M_J*M_N];
    lp = glp_create_prob();
    glp_set_prob_name(lp, M_prefix.c_str());
    glp_set_obj_dir(lp, GLP_MIN);
    glp_add_rows(lp, M_M*M_J);
    for( int m = 0; m < M_M; ++m )
    {
        for( int j = 1; j <= M_J; ++j )
        {
            glp_set_row_name(lp, m*M_J+j, (boost::format("x_%1%%2%")%m%j).str().c_str() );
            glp_set_row_bnds(lp, m*M_J+j, GLP_DB, res[m][j-1]*(1-M_tol), res[m][j-1]*(1+M_tol));
        }
    }
    glp_add_cols(lp, M_N);
    for( int n = 1; n <= M_N; ++n)
    {
        glp_set_col_name(lp, n, (boost::format("p_%1%")%n).str().c_str() );
        glp_set_col_bnds(lp, n, GLP_LO, 0.0, 0.0);
        glp_set_obj_coef(lp, n, 1.0);
    }
    int k = 1;
    for( int m = 0; m < M_M; ++m)
    {
        for( int j = 1; j <= M_J; ++j)
        {
            for( int n = 1; n <= M_N; ++n)
            {
                ia[k] = m*M_J+j;
                ja[k] = n;
                ar[k] = eval[m][j-1][n-1];
                k++;
            }
        }
    }
    glp_load_matrix(lp, M_M*M_J*M_N, ia, ja, ar);

    glp_simplex(lp, NULL);

    std::vector<int> indexes;
    for( int n = 1; n <= M_N; ++n)
    {
        auto p = glp_get_col_prim(lp, n);
        // Feel::cout << "p = " << p << std::endl;
        if( p > 1e-12 )
        {
            M_weights.push_back(p);
            indexes.push_back(n-1);
        }
    }
    glp_delete_prob(lp);

    return indexes;
}

template<typename RangeType>
double
EmpiricalQuadrature<RangeType>::evaluate( int m )
{
    tic();
    double res = 0.0, loc = 0.0;
    for( int i = 0; i < M_weights.size(); ++i )
    {
        auto const& elt = unwrap_ref( M_points[i].second );
        M_exprevals[m]->update( elt );
        // M_ctx->update( elt );
        // M_evaluator->update( vf::mapgmc( M_ctx ) );
        loc += M_weights[i]*M_exprevals[m]->eval(M_points[i].first, M_component); //M_evaluator->evalq(M_component,0,M_points[i].first)*M_ctx->J(M_points[i].first);
    }
    mpi::all_reduce(Environment::worldComm().globalComm(), loc, res, std::plus<double>());
    toc("compute");
    return res;
}

template<typename EltT, typename ExprT>
ExpressionEvaluator<EltT, ExprT>::ExpressionEvaluator( range_type const& r, expr_type const& ex )
    : super(r),
      M_expr(ex)
{}

template<typename EltT, typename ExprT>
void
ExpressionEvaluator<EltT, ExprT>::init(int o) override
{
    auto const eltForInit = boost::unwrap_ref(*boost::get<1>(M_range));
    auto q = _Q(o);

    using iim_type = vf::detail::integrate_im_type<decltype(M_range),decltype(M_expr),decltype(q),decltype(q)>;
    auto ims = iim_type::im(q,q,M_expr);
    auto im = ims.first;
    auto pts = im.points();
    auto gm = eltForInit.gm();
    auto geopc = gm->preCompute( pts );

    M_ctx = gm->template context<vm::POINT|vm::JACOBIAN|expr_type::context>( eltForInit, geopc );
    M_evaluator = std::make_shared<evaluator_type>(M_expr.evaluator( mapgmc(M_ctx) ) );
    M_weights = im.weights();
}

template<typename EltT, typename ExprT>
void
ExpressionEvaluator<EltT, ExprT>::update(element_type elt) override
{
    M_ctx->update( elt );
    M_evaluator->update( vf::mapgmc( M_ctx ) );
}

template<typename EltT, typename ExprT>
int
ExpressionEvaluator<EltT, ExprT>::nPoints() override
{
    return M_ctx->nPoints();
}

template<typename EltT, typename ExprT>
double
ExpressionEvaluator<EltT, ExprT>::weight(int i) override
{
    return M_weights[i];
}

template<typename EltT, typename ExprT>
double
ExpressionEvaluator<EltT, ExprT>::eval(int q, int comp) override
{
    return M_evaluator->evalq(comp,0,q)*M_ctx->J(q);
}

template<typename EltT, typename ExprT>
int
ExpressionEvaluator<EltT, ExprT>::order() override
{
    return M_expr.polynomialOrder();
}

#endif
