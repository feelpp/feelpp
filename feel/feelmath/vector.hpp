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

#include <cmath>
#include <vector>



namespace Feel {

//!
//! compute the exp  component wise of the vector \p v 
//!
template<typename T>
FEELPP_EXPORT std::vector<T>
exp( std::vector<T> const& v )
{
    std::vector<T> result( v.size() );
    Eigen::VectorXd v_mapped = Eigen::VectorXd::Map(&v.front(), v.size());
    Eigen::VectorXd res_mapped = Eigen::VectorXd::Map(&result.front(), result.size());
    res_mapped = v_mapped.array().exp();
    return result;
}

//!
//! compute the log (base e) component wise of the vector \p v 
//!
template<typename T>
FEELPP_EXPORT std::vector<T>
log( std::vector<T> const& v )
{
    std::vector<T> result;
    result.reserve( v.size() );
    for( auto const& i : v )
        result.push_back( std::log(i) );
    //Eigen::VectorXd v_mapped = Eigen::VectorXd::Map(&v.front(), v.size());
    //Eigen::VectorXd res_mapped = Eigen::VectorXd::Map(&result.front(), result.size());
    //res_mapped = v_mapped.array().log();
    return result;
}

//!
//! compute the log (base 10) component wise of the vector \p v 
//!
template<typename T>
FEELPP_EXPORT std::vector<T>
log10( std::vector<T> const& v )
{
    std::vector<T> result( v.size() );
    Eigen::VectorXd v_mapped = Eigen::VectorXd::Map(&v.front(), v.size());
    Eigen::VectorXd res_mapped = Eigen::VectorXd::Map(&result.front(), result.size());
    res_mapped = v_mapped.array().log10();
    return result;
}



}
