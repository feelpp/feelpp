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
//! @date 25 Nov 2017
//! @copyright 2017 Feel++ Consortium
//!
#ifndef FEELPP_FEELMODELS_HDG_OPTIONS_HPP
#define FEELPP_FEELMODELS_HDG_OPTIONS_HPP 1


#include <feel/feelmodels/modelcore/options.hpp>


namespace Feel {

namespace FeelModels {

//!
//! Mixed poisson command line options
//!
po::options_description makeMixedPoissonOptions( std::string prefix = "hdg.poisson" );
po::options_description makeMixedPoissonLibOptions( std::string prefix = "hdg.poisson" );

//!
//! options for mixed elasticity applications
//! @arg prefix name of the option section 
//! @code
//! Environment env( _argc=..., _argv=..., _desc = makeMixedElasticityOptions( <prefix>, ... );
//! @endcode
//!
po::options_description makeMixedElasticityOptions( std::string prefix = "hdg.elasticity" );

po::options_description makeMixedElasticityLibOptions( std::string prefix = "hdg.elasticity" );

//!
//! options for poroelasticity model
//!
po::options_description makeMixedPoissonElasticityOptions( std::string prefix = "hdg.poroelasticity" );

po::options_description makeMixedPoissonElasticityLibOptions( std::string prefix = "hdg.poroelasticity" );

}

}


#endif
