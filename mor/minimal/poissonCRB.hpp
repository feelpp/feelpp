//! -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*-
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
//! @author <you>
//! @date 15 Jun 2017
//! @copyright 2017 Feel++ Consortium
//!
//!

#ifndef FEELPP_CRB_MINIMAL_HPP
#define FEELPP_CRB_MINIMAL_HPP

#include <feel/feelcrb/modelcrbbase.hpp>
#include <feel/options.hpp>
#include <feel/feelmodels/modelproperties.hpp>

namespace Feel
{

po::options_description
makeOptions()
{
    FEELPP_EXPORT po::options_description options( "Poisson" );
    options.add_options()
        ( "poisson.filename", Feel::po::value<std::string>()->default_value("poisson.json"),
          "json file containing application parameters")
        ( "poisson.gamma", po::value<double>()->default_value( 1e4 ), "penalisation term" )
        ( "poisson.kappa", po::value<double>()->default_value(1), "conductivity" )
        ( "poisson.flux", po::value<double>()->default_value(1), "flux intensity at the base" )
        ( "poisson.N", po::value<int>()->default_value(100), "number of basis function to use" )
        ;
    return options;
}

FEELPP_EXPORT AboutData
makeAbout( std::string const& str = "poissoncrbmodel" )
{
    AboutData about( str.c_str() );
    return about;
}

class FEELPP_EXPORT Poisson : public ModelCrbBase<ParameterSpaceX,
                                                  FunctionSpace<Mesh<Simplex<2> >,
                                                                bases<Lagrange<1> > > >
{
public:
    using super_type = ModelCrbBase<ParameterSpaceX,
                                    FunctionSpace<Mesh<Simplex<2> >,
                                                  bases<Lagrange<1> > > >;
    using value_type = double;
    using element_type = super_type::element_type;
    using parameter_type = super_type::parameter_type;
    using mesh_type = super_type::mesh_type;
    using mesh_ptrtype = super_type::mesh_ptrtype;
    using prop_type = ModelProperties;
    using prop_ptrtype = std::shared_ptr<prop_type>;

private:
    mesh_ptrtype M_mesh;
    prop_ptrtype M_modelProps;

public:
    Poisson();
    int Qa();
    int Nl();
    int Ql( int l );

    void initModel();
    value_type
    output( int output_index, parameter_type const& mu , element_type& u, bool need_to_solve=false);
    boost::tuple<beta_vector_light_type, std::vector<beta_vector_light_type> >
    computeBetaQ( parameter_type const& mu );

}; // Poisson class

} // Feel namespace

#endif
