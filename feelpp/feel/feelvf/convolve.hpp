/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileContributor: Vincent Chabannes <vincent.chabannes@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2015-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_VF_CONVOLVE_H
#define FEELPP_VF_CONVOLVE_H

#include <feel/feelvf/expr.hpp>
#include <feel/feelvf/field.hpp>
//#include <feel/feelvf/integrator.hpp>
#include <feel/feelvf/integrate.hpp>

namespace Feel {

template <typename ... Ts>
auto convolve( Ts && ... varg )
{
    auto args = NA::make_arguments( std::forward<Ts>(varg)... );
    auto && range = args.get(_range);
    auto expr = Feel::vf::toValueExpression( args.get(_expr) );
    auto && space = args.get(_space);
    auto && quad = args.get_else( _quad, quad_order_from_expression );
    auto && quad1 = args.get_else( _quad1, quad_order_from_expression );
    GeomapStrategyType geomap = args.get_else( _geomap,GeomapStrategyType::GEOMAP_OPT );
    bool use_harts = args.get_else( _use_harts, false );
    bool verbose = args.get_else( _verbose, false );

    using _integrate_helper_type = Feel::detail::integrate_type<decltype(expr),decltype(range),decltype(quad),decltype(quad1)>;
    auto && quadptloc = args.get_else(_quadptloc, typename _integrate_helper_type::_quadptloc_ptrtype{});

    using space_type = std::decay_t<decltype(space)>;

    using T = double;
    constexpr int nDim = decay_type<space_type>::nDim;
    std::vector<Eigen::Matrix<T, nDim,1>> X,Z;
    Eigen::Matrix<T, nDim,1> Y;

    LOG(INFO) << "convolve::dof:" << space->dof()->dofPoints().size();
    size_type nLocalDofWithoutGhost = space->nLocalDofWithoutGhost();
    X.resize( nLocalDofWithoutGhost );
    for( auto const& d :  space->dof()->dofPoints() )
    {
        size_type dofId = d.first;
        if ( space->dof()->dofGlobalProcessIsGhost( dofId ) )
            continue;
        DCHECK( dofId < nLocalDofWithoutGhost ) << "dofId " << dofId << " should be less than " << nLocalDofWithoutGhost;
        auto const& y = d.second.template get<0>();
        for ( int i = 0; i < nDim; ++i )
            Y(i)=y(i);
        X[dofId] = Y;
    }
    DLOG(INFO) << "convolve::X(" << X.size() << ")=" << X;

    auto const& wc = space->worldComm();
    //mpi::all_gather( wc.localComm(), X.data(), X.size(), Z  );
    if ( wc.localSize() > 1 )
    {
        std::vector<std::vector<Eigen::Matrix<T, nDim,1>>> ZB;
        mpi::all_gather( wc.localComm(), X, ZB  );
        Z.reserve( space->nDof() );
        for ( size_type k=0;k<ZB.size();++k )
            Z.insert( Z.end(), ZB[k].begin(), ZB[k].end() );
        CHECK( Z.size() == space->nDof() ) << "Z.size=" << Z.size() << " should be "<< space->nDof();
    }
    else
    {
        Z = std::move( X );
    }
    DLOG(INFO) << "convolve::Z(" << Z.size() << ")=" << Z;
    auto ret = integrate( _range=range, _expr=expr, _quad=quad, _geomap=geomap,
                          _quad1=quad1, _verbose=verbose,_quadptloc=quadptloc).evaluate( Z, true );
    DLOG(INFO) << "convolve::ret(" << ret.size() <<")=" << ret << std::endl;
    auto v = space->element();
    size_type firstDofGlobalCluster = space->dof()->firstDofGlobalCluster();
    for( int i = 0; i < nLocalDofWithoutGhost; ++i )
        v(i) = ret[firstDofGlobalCluster+i]( 0, 0);
    sync(v);
    return v;
}

}

#endif
