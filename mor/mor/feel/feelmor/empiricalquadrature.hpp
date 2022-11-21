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

#include <feel/feelmor/parameterspace.hpp>
#include <feel/feelvf/expressionevaluator.hpp>
#include <feel/feelopt/glpk.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <numeric>

#if defined(FEELPP_HAS_GLPK_H)
#include <glpk.h>
#endif /* FEELPP_HAS_GLPK_H */

using namespace Feel;
using namespace Feel::opt;

template<typename RangeType>
class EmpiricalQuadrature
{
public:
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

    template<typename ExprT, typename FctT>
    using fct_type = typename expressionevalnl_type<ExprT, FctT>::update_function_type;

public:
    EmpiricalQuadrature( range_type const& range,
                         parameterspace_ptrtype const& Dmu,
                         std::string const& prefix = "" );

    template<typename ExprT>
    void addExpression(ExprT& ex, int comp = 0);
    template<typename ExprT>
    void addExpression(ExprT& ex, parameterelement_type& mu, int comp = 0);
    template<typename ExprT, typename FctT>
    void addExpression(ExprT& ex, parameterelement_type& mu, FctT& u, fct_type<ExprT,FctT> const& f, int comp = 0);

    void saveDatabase();
    bool loadDatabase();
    void initPoints();
    int offline();
    double evaluate( parameterelement_type const& mu, int m = 0);
    double evaluateOffline( parameterelement_type const& mu, int m = 0);

    int nbExpressions() const { return M_exprevals.size(); }
    int dimension() const { return M_weights.size(); }
    int order() const { return M_order; }
    expressionevalbase_vectype const expressionEvaluators() const { return M_exprevals; }
    parameterelement_type const sample(int i) { return M_trainset->at(i); }
    void setDbFilename( std::string const& name ) { M_dbFilename = name; }

private:
    void solveLP( evaluation_type const& eval, local_result_type const& res );

private:
    range_type M_range;
    expressionevalbase_vectype M_exprevals;
    std::vector<std::vector<double> > M_evaluations;
    std::vector<double> M_weights;
    std::vector<std::pair<int,element_type> >  M_points;
    std::vector<int> M_indexes;
    parameterspace_ptrtype M_Dmu;
    sampling_ptrtype M_trainset;

    std::string M_prefix;
    int M_verbosity;
    bool M_dbLoad;
    std::string M_dbFilename;
    double M_tol;
    int M_M;
    int M_J;
    int M_numElts;
    int M_N;
    int M_order;
    int M_maxOrder;
    double M_tolZero;
};

template<typename RangeType>
EmpiricalQuadrature<RangeType>::EmpiricalQuadrature( range_type const& range,
                                                     parameterspace_ptrtype const& Dmu,
                                                     std::string const& prefix )
    :
    M_range(range),
    M_Dmu(Dmu),
    M_prefix(prefix),
    M_verbosity(ioption(_prefix=M_prefix,_name="eq.verbosity")),
    M_dbLoad(boption(_prefix=M_prefix,_name="eq.db.load")),
    M_dbFilename(soption(_prefix=M_prefix,_name="eq.db.filename")),
    M_tol(doption(_prefix=M_prefix,_name="eq.tolerance")),
    M_M(0),
    M_J(ioption(_prefix=M_prefix,_name="eq.sampling-size")),
    M_numElts(0),
    M_N(0),
    M_order(ioption(_prefix=M_prefix, _name="eq.order")),
    M_maxOrder(ioption(_prefix=M_prefix, _name="eq.max-order")),
    M_tolZero(doption(_prefix=M_prefix, _name="eq.tolerance-zero"))
{
}

template<typename RangeType>
template<typename ExprT>
void
EmpiricalQuadrature<RangeType>::addExpression( ExprT& ex, int comp )
{
    auto exev = std::make_shared<expressioneval_type<ExprT>>(M_range, ex, comp);
    M_exprevals.push_back(exev);
    M_M++;
}

template<typename RangeType>
template<typename ExprT>
void
EmpiricalQuadrature<RangeType>::addExpression( ExprT& ex, parameterelement_type& mu, int comp )
{
    auto exev = std::make_shared<expressionevalparam_type<ExprT>>(M_range, ex, mu, comp);
    M_exprevals.push_back(exev);
    M_M++;
}

template<typename RangeType>
template<typename ExprT, typename FctT>
void
EmpiricalQuadrature<RangeType>::addExpression( ExprT& ex, parameterelement_type& mu,
                                               FctT& u, fct_type<ExprT,FctT> const& f, int comp )
{
    auto exev = std::make_shared<expressionevalnl_type<ExprT, FctT>>(M_range, ex, mu, u, comp);
    exev->M_fct = f;
    M_exprevals.push_back(exev);
    M_M++;
}

template<typename RangeType>
void
EmpiricalQuadrature<RangeType>::saveDatabase()
{
    fs::path dbFilenamePath = M_dbFilename;
    if( Environment::isMasterRank() )
    {
        if( fs::exists(dbFilenamePath) )
        {
            if( !fs::is_directory(dbFilenamePath) )
                throw std::exception();
        }
        else
            fs::create_directory(dbFilenamePath);
    }
    Environment::worldComm().barrier();
    dbFilenamePath /= "db_" + std::to_string(Environment::rank());
    std::ofstream os(dbFilenamePath.string(), std::ios::binary);
    boost::archive::binary_oarchive oar(os);
    oar << M_indexes;
    oar << M_weights;
}

template<typename RangeType>
bool
EmpiricalQuadrature<RangeType>::loadDatabase()
{
    fs::path dbFilenamePath = M_dbFilename;
    dbFilenamePath /= "db_" + std::to_string(Environment::rank());
    if( fs::exists(dbFilenamePath) )
    {
        std::ifstream is(dbFilenamePath.string(), std::ios::binary);
        boost::archive::binary_iarchive iar(is);
        iar >> M_indexes;
        iar >> M_weights;
        return true;
    }
    else
        return false;
}

template<typename RangeType>
void
EmpiricalQuadrature<RangeType>::initPoints()
{
    std::cout << "[" << Environment::rank() << "] non zero : " << M_indexes.size() << std::endl;
    int n = 0;
    auto it = M_indexes.begin();
    auto itEnd = M_indexes.end();
    int nPts = M_exprevals[0]->nPoints();
    for( auto const& eltWrap : M_range )
    {
        auto const& elt = unwrap_ref( eltWrap );
        if ( elt.processId() != Environment::rank() )
            continue;
        for ( uint16_type q = 0; q < nPts; ++q )
        {
            if( it != itEnd && n == *it )
            {
                ++it;
                M_points.push_back(std::make_pair(q,eltWrap));
            }
            n++;
        }
    }
}

template<typename RangeType>
int
EmpiricalQuadrature<RangeType>::offline()
{
    if( M_exprevals.size() == 0 )
        return 0;

    if( M_order < 0 )
    {
        int max_order = std::accumulate(M_exprevals.begin(), M_exprevals.end(), 0,
                                        [](auto const& a, auto const& b)
                                            {
                                                return std::max(a,b->order());
                                            });
        M_order = std::min(max_order,M_maxOrder);
    }
    for( auto const& ee : M_exprevals )
        ee->init(M_order);

    if( M_dbLoad && this->loadDatabase() )
    {
        this->initPoints();
        return 0;
    }

    M_numElts = nelements(M_range);
    if( M_numElts <= 0 )
        return 0;

    M_trainset = M_Dmu->sampling();
    M_trainset->randomize(M_J);

    tic();
    auto const eltForInit = boost::unwrap_ref(*boost::get<1>(M_range));

    int nPts = M_exprevals[0]->nPoints();
    M_N = M_numElts*nPts;
    evaluation_type eval(M_M, std::vector<std::vector<double> >(M_J, std::vector<double>(M_N, 0.0)));
    local_result_type res(M_M, std::vector<double>(M_J, 0.0) );
    std::vector<int> muToRemove;
    int m = 0, j = 0, n = 0;
    Feel::cout << "offline with M = " << M_M << " J = " << M_J << " M_N = " << M_N
               << " (elts=" << M_numElts << " nPoints=" << nPts << ")"
               << " order = " << M_order << std::endl;
    for( j = 0; j < M_J; ++j )
    {
        bool b = true;
        for( m = 0; m < M_M; ++m )
            b = b && M_exprevals[m]->update(M_trainset->at(j));
        if( !b )
        {
            muToRemove.push_back(j);
            continue;
        }
        n = 0;
        for( auto const& eltWrap : M_range )
        {
            auto const& elt = unwrap_ref( eltWrap );
            if ( elt.processId() != Environment::rank() )
                continue;

            for( m = 0; m < M_M; ++m )
            {
                M_exprevals[m]->update( eltWrap );
                for ( uint16_type q = 0; q < nPts; ++q )
                {
                    eval[m][j][n+q] = M_exprevals[m]->eval(q, M_exprevals[m]->component());
                    res[m][j] += M_exprevals[m]->weight(q)*eval[m][j][n+q];
                }
            }
            n += nPts;
        }
    }
    for( auto it = muToRemove.rbegin(); it != muToRemove.rend(); ++it )
    {
        M_trainset->erase(std::next(M_trainset->begin(),*it));
        for( m = 0; m < M_M; ++m )
        {
            eval[m].erase(std::next(eval[m].begin(),*it));
            res[m].erase(std::next(res[m].begin(),*it));
        }
    }
    M_J -= muToRemove.size();
    toc("evaluations");

    // for( int m = 0; m < M_M; ++m )
    // {
    //     double b = 0.0;
    //     mpi::all_reduce(Environment::worldComm().globalComm(), res[m][0], b, std::plus<double>());
    //     Feel::cout << std::setprecision(10) << "b[" << m << "] = " << b << std::endl;
    //     std::cout << "[" << Environment::rank() << "] res = " << res[m][0] << std::endl;
    // }

    try{
        solveLP(eval,res);
    }
    catch( FeelGlpkException& e )
    {
        Feel::cout << std::string(e.what()) << std::endl;
        return e.M_e;
    }

    this->initPoints();
    this->saveDatabase();

    return 0;
}

template<typename RangeType>
void
EmpiricalQuadrature<RangeType>::solveLP( evaluation_type const& eval, local_result_type const& res )
{
    auto glpk = OptimizationLinearProgramming( GLP_MIN, M_prefix );

    tic();
    for( int m = 1; m <= M_M; ++m )
        for( int j = 1; j <= M_J; ++j )
            glpk.addRow((boost::format("x_%1%%2%")%m%j).str(), GLP_DB, res[m-1][j-1]-M_tol, res[m-1][j-1]+M_tol );

    for( int n = 1; n <= M_N; ++n)
        glpk.addColumn( (boost::format("p_%1%")%n).str(), 1.0, GLP_LO, 0.0, 0.0 );

    std::vector<std::vector<double> > eval2(M_M*M_J);
    for(int m = 0; m < M_M; ++m)
        for(int j = 0; j < M_J; ++j)
            for(int n = 0; n < M_N; ++n)
                eval2[m*M_J+j].push_back(eval[m][j][n]);
    glpk.setMatrix(eval2);
    toc("init lp");

    tic();
    int e = glpk.solve();
    if( e != 0 )
    {
        Feel::cout << tc::red << glpk.check(e) << tc::reset << std::endl;
        throw FeelGlpkException(e);;
    }
    toc("solve lop");

    M_indexes.clear();
    for( int n = 1; n <= M_N; ++n)
    {
        double p = glpk.getColumnPrimalValue(n);
        if( p > M_tolZero )
        {
            M_weights.push_back(p);
            M_indexes.push_back(n-1);
        }
    }
}

template<typename RangeType>
double
EmpiricalQuadrature<RangeType>::evaluate( parameterelement_type const& mu, int m )
{
    tic();
    double res = 0.0, loc = 0.0;
    M_exprevals[m]->update(mu);
    toc("update", M_verbosity > 1 );
    tic();
    for( int i = 0; i < M_weights.size(); ++i )
    {
        M_exprevals[m]->update( M_points[i].second );
        loc += M_weights[i]*M_exprevals[m]->eval(M_points[i].first, M_exprevals[m]->component());
    }
    mpi::all_reduce(Environment::worldComm().globalComm(), loc, res, std::plus<double>());
    toc("compute", M_verbosity > 0 );
    return res;
}

template<typename RangeType>
double
EmpiricalQuadrature<RangeType>::evaluateOffline( parameterelement_type const& mu, int m )
{
    double res = 0.0, loc = 0.0;
    tic();
    for( auto const& eltWrap : M_range )
    {
        auto const& elt = unwrap_ref( eltWrap );
        if ( elt.processId() != Environment::rank() )
            continue;

        M_exprevals[m]->update( eltWrap );

        for ( uint16_type q = 0; q < M_exprevals[m]->nPoints(); ++q )
            loc += M_exprevals[m]->weight(q)*M_exprevals[m]->eval(q, M_exprevals[m]->component());
    }

    mpi::all_reduce(Environment::worldComm().globalComm(), loc, res, std::plus<double>());

    toc("evaluateOffline");
    return res;
}

#endif
