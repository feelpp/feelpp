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
//! @date 10 Jun 2017
//! @copyright 2017 Feel++ Consortium
//!

#ifndef FEELPP_CRBDATA_HPP
#define FEELPP_CRBDATA_HPP 1

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <Eigen/Core>
#include <feel/feelmor/parameterspace.hpp>



namespace Feel {

using vectorN_type=Eigen::VectorXd;
using matrixN_type=Eigen::MatrixXd;
namespace RB {
typedef boost::tuple< std::vector<vectorN_type> , std::vector<vectorN_type> , std::vector<vectorN_type>, std::vector<vectorN_type> > solutions_tuple;
typedef boost::tuple< std::vector<double>,double,double , std::vector< std::vector< double > > , std::vector< std::vector< double > > > upper_bounds_tuple;
typedef boost::tuple< double,double > matrix_info_tuple; //conditioning, determinant
}

class CRBSolutions : public RB::solutions_tuple
{
public:
    
};
class CRBMatrixInfo : public RB::matrix_info_tuple
{
public:
    double conditioning() const { return boost::get<0>( *this ); }
    double determinant() const { return boost::get<1>( *this ); }
};
class CRBUpperBounds : public RB::upper_bounds_tuple
{
public:
    
};

//!
//! data structure holding the results of a RB online code
//!
class CRBResults
    : public boost::tuple<std::vector<double>,double, RB::solutions_tuple, RB::matrix_info_tuple, double, double, RB::upper_bounds_tuple > 
{
    using super = boost::tuple<std::vector<double>,double, RB::solutions_tuple, RB::matrix_info_tuple, double, double, RB::upper_bounds_tuple >;
public:
    CRBResults() = default;
    CRBResults( CRBResults const& ) = default;
    CRBResults( CRBResults && ) = default;
    CRBResults& operator=( CRBResults const& ) = default;
    CRBResults& operator=( CRBResults && ) = default;
    using parameter_type = ParameterSpaceX::Element;
    void setParameter( parameter_type const& mu )
        {
            M_mu = mu;
        }
    parameter_type const& parameter() const { return M_mu; }
    //!
    //! @return the output
    //!
    double output() const { return boost::get<0>( *this ).back(); }
    double errorbound() const  { return boost::get<0>( boost::get<6>( *this ) ).back(); }
    vectorN_type const& coefficients() const { return boost::get<0>(boost::get<2>( *this )).back(); }

    CRBResults ( super const& s )
        :
        super( s )
        {}
    CRBResults ( super && s )
        :
        super( std::move(s) )
        {}
private:
    parameter_type M_mu;
};


}

#endif
