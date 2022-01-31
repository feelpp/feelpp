/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

 This file is part of the Feel library

 Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 Date: 2009-11-13

 Copyright (C) 2009-2014 Feel++ Consortium

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
/**
 \file heat2d.hpp
 \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 \date 2009-11-13
 */
#ifndef FEELPP_HEAT2D_HPP
#define FEELPP_HEAT2D_HPP 1

#include <boost/timer.hpp>

#include <feel/options.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/gmsh.hpp>

#include <feel/feelmor/modelcrbbase.hpp>


namespace Feel
{

FEELPP_EXPORT po::options_description
makeHeat2DOptions()
{
    po::options_description heat2doptions( "Heat2D options" );
    heat2doptions.add_options()
        ( "gamma", po::value<double>()->default_value( 1e4 ), "penalisation term" );
        //( "k1", po::value<double>()->default_value( 0.1 ), "k1" )
        //( "k2", po::value<double>()->default_value( 0.1 ), "k2" )
    return heat2doptions;
}
FEELPP_EXPORT AboutData
makeHeat2DAbout( std::string const& str = "heat2d" )
{
    Feel::AboutData about( /*AppName  */ str.c_str(),
                           /*ProgName */ str.c_str(),
                           /*Version  */ "0.1",
                           /*ShortDesc*/ "2D Heat Benchmark",
                           /*Licence  */ Feel::AboutData::License_GPL,
                           /*Copyright*/ "Copyright (c) 2009-2014 Feel++ Consortium" );

    about.addAuthor( "Vincent HUBER", "developer", "vincent.huber@cemosis.fr", "" );
    return about;
}

/**
 * \class Heat2D
 * \brief brief description
 *
 * @author Vincent HUBER
 * @see
 */
class FEELPP_EXPORT  Heat2D : public ModelCrbBase<ParameterSpaceX, decltype(Pch<3>(Mesh<Simplex<2>>::New()))>
{
public:
    using super = ModelCrbBase<ParameterSpaceX, decltype(Pch<3>(Mesh<Simplex<2>>::New()))>;

    Heat2D() : super( "heat2d" ) {}

    //! initialisation of the model
    void initModel();

    //! compute beta coefficients
    super::betaq_type computeBetaQ( parameter_type const& mu );

    /**
     * Given the output index \p output_index and the parameter \p mu, return
     * the value of the corresponding FEM output
     */
    value_type output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve=false);
};


}

#endif /* FEELPP_HEAT2D_HPP */
