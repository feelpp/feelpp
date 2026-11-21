/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 *
 * SPDX-FileCopyrightText: 2010 Universite Joseph Fourier (Grenoble I)
 * SPDX-FileCopyrightText: 2011-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
/**
   \file norm.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2010-04-17
 */
#ifndef __norm_H
#define __norm_H 1

namespace Feel
{
namespace vf
{
/**
 * inner product of the left and right expressions
 *
 * \note it handles vectorial expressions
 *
 * \return the inner product \f$(l * r^T)\f$
 */
template<typename ExprL, typename ExprR>
inline
auto
outer_prod( ExprL l, ExprR r ) -> decltype( l * trans( r ) )
{
    return l * trans( r );
}

/**
 * double dot product of the left and right expressions
 *
 * \note it handles matricial(rank 2 tensors) expressions
 *
 * \return the double dot product \f$(l : r)\f$
 */
template<typename ExprL, typename ExprR>
requires requires
{
    typename ExprL::evaluate_type;
    typename ExprR::evaluate_type;
} &&
( ExprL::evaluate_type::RowsAtCompileTime == ExprR::evaluate_type::RowsAtCompileTime ) &&
( ExprL::evaluate_type::ColsAtCompileTime == ExprR::evaluate_type::ColsAtCompileTime )
inline
auto
ddot( ExprL l, ExprR r ) -> decltype( trace( trans( l ) * r ) )
{
    return trace( trans( l ) * r );
}


} // vf


} // feel
#endif /* __norm_H */
