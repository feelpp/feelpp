//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
//!
//! Author(s) : 
//!     Thibaut Metivet <thibaut.metivet@inria.fr>
//!
//! This library is free software; you can redistribute it and/or
//! modify it under the terms of the GNU Lesser General Public
//! License as published by the Free Software Foundation; either
//! version 2.1 of the License, or (at your option) any later version.
//!
//! This library is distributed in the hope that it will be useful,
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//! Lesser General Public License for more details.
//!
//! You should have received a copy of the GNU Lesser General Public
//! License along with this library; if not, write to the Free Software
//! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//!
//! @file simplexeikonalsolver.hpp
//! @author Thibaut Metivet <thibaut.metivet@inria.fr>
//! @date 12 Dec 2019
//! @copyright 2019 Feel++ Consortium
//! @copyright 2019 INRIA
//!

#ifndef _SIMPLEX_EIKONAL_SOLVER_HPP
#define _SIMPLEX_EIKONAL_SOLVER_HPP 1

#include <Eigen/Core>
#include <Eigen/QR>
#include "eigenmap.hpp"
#include "fastmarchingdofstatus.hpp"

namespace Feel {

template< typename FunctionSpaceType >
class SimplexEikonalSolver
{
    public:
        //--------------------------------------------------------------------//
        // Functionspace and mesh
        typedef FunctionSpaceType functionspace_type;
        typedef std::shared_ptr< functionspace_type > functionspace_ptrtype;
        static const uint16_type nDofPerElt = functionspace_type::fe_type::nDof;
        typedef typename functionspace_type::element_type element_type;
        typedef typename functionspace_type::element_ptrtype element_ptrtype;
        //--------------------------------------------------------------------//
        static constexpr uint16_type nRealDim = functionspace_type::nRealDim;
        typedef typename functionspace_type::value_type value_type;
        typedef typename node<value_type>::type node_type;
        typedef typename matrix_node<value_type>::type matrix_node_type;
        //--------------------------------------------------------------------//
        typedef std::pair< size_type, value_type > pair_dof_value_type;

    public:
        SimplexEikonalSolver( functionspace_ptrtype const& space ) : M_space( space ) {}

        std::vector< pair_dof_value_type > solveEikonal( std::vector< size_type > dofCloseIds, std::vector< size_type > const& dofDoneIds, size_type eltId, element_type const& sol );

    private:
        functionspace_ptrtype M_space;
};

template< typename FunctionSpaceType >
std::vector< typename SimplexEikonalSolver< FunctionSpaceType >::pair_dof_value_type >
SimplexEikonalSolver< FunctionSpaceType >::solveEikonal( std::vector< size_type > dofCloseIds, std::vector< size_type > const& dofDoneIds, size_type eltId, element_type const& sol )
{
    auto const& dofTable = sol.functionSpace()->dof();

    std::vector< std::pair< size_type, value_type > > closeDofsValues;

    size_type nDofClose = dofCloseIds.size();
    if( nDofClose == 0 )
        return closeDofsValues;
    size_type nDofDone = dofDoneIds.size();

    auto const& pt0 = boost::get<0>( dofTable->dofPoint( dofDoneIds[0] ) );
    auto const x0 = eigenMap<nRealDim>( pt0 );
    value_type const phi0 = std::abs( sol(dofDoneIds[0]) );

    if( nDofDone == 1 )
    {
        for( size_type dofCloseId: dofCloseIds )
        {
            auto const ptClose = boost::get<0>( dofTable->dofPoint( dofCloseId ) );
            auto const xClose = eigenMap<nRealDim>( ptClose );
            closeDofsValues.emplace_back( dofCloseId, phi0 + (xClose - x0).norm() );
        }
        return closeDofsValues;
    }

    // Compute deltaPhi(i) = phi(i+1)-phi(0);
    Eigen::Matrix< value_type, Eigen::Dynamic, 1 > deltaPhi( nDofDone - 1 );
    for( size_type i = 0; i < nDofDone-1; ++i )
    {
        deltaPhi(i) = std::abs( sol(dofDoneIds[i+1]) ) - phi0;
    }

    // Compute orthonormal basis from [ x1-x0, ..., xi-x0 ] (using QR Householder)
    Eigen::Matrix< value_type, nRealDim, Eigen::Dynamic > deltaX( nRealDim, nDofDone-1 );
    for( size_type i = 0; i < nDofDone-1; ++i )
    {
        auto const& pti = boost::get<0>( dofTable->dofPoint( dofDoneIds[i+1] ) );
        auto const xi = eigenMap<nRealDim>( pti );
        deltaX.col(i) = xi - x0;
    }

    auto const deltaX_QR = deltaX.householderQr();
    // Ensures positive diagonal for R
    auto const deltaX_R_raw = deltaX_QR.matrixQR().topLeftCorner( nDofDone-1, nDofDone-1 );
    Eigen::Matrix< value_type, Eigen::Dynamic, Eigen::Dynamic > deltaX_Q_raw = deltaX_QR.householderQ();
    Eigen::DiagonalMatrix< value_type, Eigen::Dynamic > signDiagR( nDofDone - 1 );
    signDiagR.diagonal() = deltaX_R_raw.diagonal().array().sign();
    //auto const signDiagR = deltaX_R_raw.diagonal().array().sign().matrix().asDiagonal();
    Eigen::Matrix< value_type, Eigen::Dynamic, Eigen::Dynamic > deltaX_R = signDiagR * deltaX_R_raw;
    Eigen::Matrix< value_type, Eigen::Dynamic, Eigen::Dynamic > deltaX_Q = deltaX_Q_raw.leftCols( nDofDone-1 ) * signDiagR;
    
    // Compute first (nDofDone-1) gradient coefficients
    Eigen::Matrix< value_type, Eigen::Dynamic, 1 > gradCoeff( nDofDone, 1 );
    gradCoeff.topRows( nDofDone-1 ) = deltaX_R
        .transpose()
        .template triangularView<Eigen::Lower>()
        .solve( deltaPhi );
    // Compute the remaining gradient coefficient (in direction xClose)
    value_type restTo1 = 1. - gradCoeff.topRows( nDofDone - 1 ).squaredNorm();
    if( restTo1 < 0. )
    {
        //gradCoeff( nDofDone-1 ) = 0.;
        return closeDofsValues;
    }
    else
    {
        gradCoeff( nDofDone-1 ) = std::sqrt( restTo1 );
    }

    // Then try to solve for the close dof values
    for( auto dofCloseId: dofCloseIds )
    {
        auto const ptClose = boost::get<0>( dofTable->dofPoint( dofCloseId ) );
        auto const xClose = eigenMap<nRealDim>( ptClose );
        auto const deltaXClose = xClose - x0;

        Eigen::Matrix< value_type, Eigen::Dynamic, 1 > deltaXCloseCoeff( nDofDone, 1 );
        deltaXCloseCoeff.topRows( nDofDone - 1 ) = deltaX_Q.transpose() * deltaXClose;
        deltaXCloseCoeff( nDofDone - 1 ) = ( deltaXClose - ( deltaX_Q * ( deltaX_Q.transpose() * deltaXClose ).asDiagonal() ).rowwise().sum() ).normalized().dot( deltaXClose );

        // Check upwind criterion: solve deltaX^T lambda = xf = xpt-x0 - ((xpt-x0).en)/(g.en) g and check lambda_i in [0,1]
        Eigen::Matrix< value_type, Eigen::Dynamic, 1 > lambda( nDofDone-1, 1 );
        auto const xfCoeff = deltaXCloseCoeff - deltaXCloseCoeff( nDofDone - 1 ) / gradCoeff( nDofDone-1 ) * gradCoeff;
        lambda = deltaX_R
            .template triangularView<Eigen::Upper>()
            .solve( xfCoeff.topRows( nDofDone - 1 ) );
        bool upwind = ( 0. <= lambda.array() ).all() && ( lambda.array() <= 1. ).all();

        // Return solution = phi0 + grad . (xpt-x0)
        value_type val;
        if( upwind )
        {
            val =  phi0 + gradCoeff.transpose() * deltaXCloseCoeff;
        }
        else
        {
            // Compute distance to new DONE dof (dofDoneId)
            val = phi0 + (xClose - x0).norm();
        }

        closeDofsValues.emplace_back( dofCloseId, val );
    }

    return closeDofsValues;
}

} // namespace Feel

#endif // _SIMPLEX_EIKONAL_SOLVER_HPP

