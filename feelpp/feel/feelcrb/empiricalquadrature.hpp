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
#include <feel/feelvf/expressionevaluator.hpp>
#include <numeric>

#if defined(FEELPP_HAS_GLPK_H)
#include <glpk.h>
#endif /* FEELPP_HAS_GLPK_H */

using namespace Feel;

template<typename RangeType>
class EmpiricalQuadrature
{
    using range_type = RangeType;
    using element_type = typename boost::tuples::element<1_c,range_type>::type::value_type;
    using gm_type = typename element_type::type::gm_type;

    using expressionevalbase_type = ExpressionEvaluatorBase<range_type>;
    using expressionevalbase_ptrtype = std::shared_ptr<expressionevalbase_type>;
    using expressionevalbase_vectype = std::vector<expressionevalbase_ptrtype>;
    template<typename ExprT>
    using expressioneval_type = ExpressionEvaluator<range_type, ExprT>;
    template<typename ExprT>
    using expressionevalparam_type = ExpressionEvaluatorParam<range_type, ExprT>;
    template<typename ExprT, typename FctT>
    using expressionevalnl_type = ExpressionEvaluatorNonLinear<range_type, ExprT, FctT>;

    using parameterspace_type = ParameterSpace<>;
    using parameterspace_ptrtype = std::shared_ptr<parameterspace_type>;
    using parameterelement_type = parameterspace_type::element_type;
    using sampling_type = typename parameterspace_type::sampling_type;
    using sampling_ptrtype = typename parameterspace_type::sampling_ptrtype;

    using evaluation_type = std::vector<std::vector<std::vector<double> > >;
    using local_result_type = std::vector<std::vector<double> >;

    template<typename FctT>
    using fct_type = std::function<void(parameterelement_type const&, FctT&)>;

public:
    EmpiricalQuadrature( range_type const& range,
                         parameterspace_ptrtype const& Dmu,
                         int comp = 0, std::string const& prefix = "" );

    template<typename ExprT>
    void addExpression(ExprT& ex);
    template<typename ExprT>
    void addExpression(ExprT& ex, parameterelement_type& mu);
    template<typename ExprT, typename FctT>
    void addExpression(ExprT& ex, parameterelement_type& mu, FctT& u, fct_type<FctT> const& f);

    int offline();
    double evaluate( parameterelement_type const& mu, int m = 0);

    int dimension() const { return M_weights.size(); }
    expressionevalbase_vectype const expressionEvaluators() const { return M_exprevals; }

private:
    std::vector<int> solveLP( evaluation_type const& eval, local_result_type const& res );

private:
    range_type M_range;
    expressionevalbase_vectype M_exprevals;
    int M_component;
    std::vector<std::vector<double> > M_evaluations;
    std::vector<double> M_weights;
    std::vector<std::pair<int,element_type> >  M_points;
    parameterspace_ptrtype M_Dmu;
    sampling_ptrtype M_trainset;

    std::string M_prefix;
    double M_tol;
    int M_M;
    int M_J;
    int M_numElts;
    int M_N;
};

template<typename RangeType>
EmpiricalQuadrature<RangeType>::EmpiricalQuadrature( range_type const& range,
                                                     parameterspace_ptrtype const& Dmu,
                                                     int comp, std::string const& prefix )
    :
    M_range(range),
    M_component(comp),
    M_Dmu(Dmu),
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
EmpiricalQuadrature<RangeType>::addExpression( ExprT& ex)
{
    auto exev = std::make_shared<expressioneval_type<ExprT>>(M_range, ex);
    M_exprevals.push_back(exev);
    M_M++;
}

template<typename RangeType>
template<typename ExprT>
void
EmpiricalQuadrature<RangeType>::addExpression( ExprT& ex, parameterelement_type& mu)
{
    auto exev = std::make_shared<expressionevalparam_type<ExprT>>(M_range, ex, mu);
    M_exprevals.push_back(exev);
    M_M++;
}

template<typename RangeType>
template<typename ExprT, typename FctT>
void
EmpiricalQuadrature<RangeType>::addExpression( ExprT& ex, parameterelement_type& mu,
                                               FctT& u, fct_type<FctT> const& f )
{
    auto exev = std::make_shared<expressionevalnl_type<ExprT, FctT>>(M_range, ex, mu, u);
    exev->M_fct = f;
    M_exprevals.push_back(exev);
    M_M++;
}

template<typename RangeType>
int
EmpiricalQuadrature<RangeType>::offline()
{
    if( M_numElts <= 0 )
        return 1;

    tic();
    auto const eltForInit = boost::unwrap_ref(*boost::get<1>(M_range));

    int max_order = std::accumulate(M_exprevals.begin(), M_exprevals.end(), 0,
                                    [](auto const& a, auto const& b)
                                        {
                                            return std::max(a,b->order());
                                        });
    for( auto const& ee : M_exprevals )
        ee->init(max_order);

    int nPts = M_exprevals[0]->nPoints();
    M_N = M_numElts*nPts;
    evaluation_type eval(M_M, std::vector<std::vector<double> >(M_J, std::vector<double>(M_N, 0.0)));
    local_result_type res(M_M, std::vector<double>(M_N, 0.0) );
    int m = 0, j = 0, n = 0;
    for( m = 0; m < M_M; ++m )
    {
        for( j = 0; j < M_J; ++j )
        {
            M_exprevals[m]->update(M_trainset->at(j));
            n = 0;
            for( auto const& eltWrap : M_range )
            {
                auto const& elt = unwrap_ref( eltWrap );
                if ( elt.processId() != Environment::rank() )
                    continue;

                for( auto const& ee : M_exprevals )
                    ee->update( eltWrap );

                for ( uint16_type q = 0; q < nPts; ++q, ++n )
                {
                    eval[m][j][n] = M_exprevals[m]->eval(q, M_component);
                    res[m][j] += M_exprevals[m]->weight(q)*eval[m][j][n];
                }
            }
        }
    }
    toc("evaluations");

    // for( int m = 0; m < M_M; ++m )
    // {
    //     double b = 0.0;
    //     mpi::all_reduce(Environment::worldComm().globalComm(), res[m][0], b, std::plus<double>());
    //     Feel::cout << "b[" << m << "] = " << b << std::endl;
    // }

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
    toc("init lp");

    tic();
    glp_simplex(lp, NULL);
    toc("solve lop");

    std::vector<int> indexes;
    for( int n = 1; n <= M_N; ++n)
    {
        auto p = glp_get_col_prim(lp, n);
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
EmpiricalQuadrature<RangeType>::evaluate( parameterelement_type const& mu, int m )
{
    tic();
    double res = 0.0, loc = 0.0;
    for( int i = 0; i < M_weights.size(); ++i )
    {
        M_exprevals[m]->update(mu);
        M_exprevals[m]->update( M_points[i].second );
        loc += M_weights[i]*M_exprevals[m]->eval(M_points[i].first, M_component);
    }
    mpi::all_reduce(Environment::worldComm().globalComm(), loc, res, std::plus<double>());
    toc("compute");
    return res;
}

#endif
