//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
//!
//! This file is part of the Feel++ library
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
//! @file
//! @author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
//! @date 24 Jul 2017
//! @copyright 2017 Feel++ Consortium
//!

#include <Eigen/QR>
#include <vector>

namespace Feel {

//!
//! compute least square (or interpolation) fitting polynomial
//! @param xv x parameter values
//! @param yv y parameter values
//! @param order order of the least square (or interpolation) fitting polynomial
//!
//! Note:
//!   . xv.size() == order+1 then it is interpolation,
//!   . xv.size() > order+1 then it is least square fitting
//!
template<typename T>
FEELPP_EXPORT std::vector<T>
polyfit(const std::vector<T> &xv, const std::vector<T> &yv, int order = 1)
{
	Eigen::MatrixXd A(xv.size(), order+1);
    Eigen::VectorXd xv_mapped = Eigen::VectorXd::Map(&xv.front(), xv.size());
	Eigen::VectorXd yv_mapped = Eigen::VectorXd::Map(&yv.front(), yv.size());
	Eigen::VectorXd result;
    
	if (xv.size() != yv.size()) throw std::logic_error("invalid input vectors size");
	if (xv.size() < order+1) throw  std::logic_error(std::string("input vectors size must be at least >= ") + std::to_string(order+1) ) ;

	// create VDM matrix
	for (size_t i = 0; i < xv.size(); i++)
        for (size_t j = 0; j < order+1; j++)
            A(i, j) = std::pow(xv.at(i), j);
    
	// solve for linear least squares fit
	result = A.householderQr().solve(yv_mapped);
    
    std::vector<T> coeff( order+1 );
	for (size_t i = 0; i < order+1; i++)
		coeff[i] = result[i];
    
    return coeff;
}



} // Feel
