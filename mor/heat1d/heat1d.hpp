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
//! @date 14 Jun 2017
//! @copyright 2017 Feel++ Consortium
//!
#ifndef FEELPP_HEAT1D_HPP
#define FEELPP_HEAT1D_HPP 1

#include <feel/options.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/gmsh.hpp>

#include <feel/feelcrb/modelcrbbase.hpp>


namespace Feel
{

FEELPP_EXPORT po::options_description
makeHeat1DOptions()
{
    po::options_description heat1doptions( "Heat1D options" );
    // heat1doptions.add_options()
    // ( "mu1", po::value<double>()->default_value( 0.2 ), "mu1" )
    // ( "mu2", po::value<double>()->default_value( 0.2 ), "mu2" )
    // ( "mu3", po::value<double>()->default_value(-1.0 ), "mu3" )
    // ( "mu4", po::value<double>()->default_value( 0.1 ), "mu4" )
    // ;
    return heat1doptions;
}
FEELPP_EXPORT AboutData
makeHeat1DAbout( std::string const& str = "heat1d" )
{
    Feel::AboutData about( /*AppName  */ str.c_str(),
                           /*ProgName */ str.c_str(),
                           /*Version  */ "0.1",
                           /*ShortDesc*/ "1D Heat Benchmark",
                           /*Licence  */ Feel::AboutData::License_GPL,
                           /*Copyright*/ "Copyright (c) 2009-2014 Feel++ Consortium" );

    about.addAuthor( "Christophe Prud'homme", "developer", "christophe.prudhomme@feelpp.org", "" );
    return about;
}

/**
 * \class Heat1D
 * \brief brief description
 *
 * @author Christophe Prud'homme
 * @see
 */
class FEELPP_EXPORT Heat1D : public ModelCrbBase<ParameterSpaceX, decltype(Pch<5>(Mesh<Simplex<1>>::New()))>
{
public:

    typedef ModelCrbBase<ParameterSpaceX, decltype(Pch<5>(Mesh<Simplex<1>>::New()))> super_type;

    Heat1D() : super_type( "heat1d" ) {}

    void initBetaQ();
        
    super_type::betaq_type computeBetaQ( parameter_type const& mu );

    //! initialisation of the model
    void initModel();

    beta_vector_light_type beta;
    /**
     * Given the output index \p output_index and the parameter \p mu, return
     * the value of the corresponding FEM output
     */
    value_type output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve=false);
};



}

#endif /* FEELPP_HEAT1D_HPP */
