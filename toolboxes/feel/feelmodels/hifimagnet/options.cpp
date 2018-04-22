/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

   This file is part of the Feel library

   Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   Author(s): Cecile Daversin <cecile.daversin@lncmi.cnrs.fr>
   Date: 2011-16-12

   Copyright (C) 2008-2010 Universite Joseph Fourier (Grenoble I)
   Copyright (C) CNRS

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
   \file options.cpp
   \author Cecile Daversin <cecile.daversin@lncmi.cnrs.fr>
   \date 2013-21-05
*/

#include<feel/options.hpp>
#include<feel/feelmodels/hifimagnet/options.hpp>

Feel::po::options_description
HifiMagnetOptions()
{
    Feel::po::options_description hifimagnet_options("hifimagnet general options");
    hifimagnet_options.add_options()
        // Dimension
        ("dim", Feel::po::value<int>()->default_value( 2 ), "dimension")
        // Geofile are meshfile inputs
        ("geofile", Feel::po::value<std::string>()->default_value( "" ), "name of the geofile input (geo or mesh)")
        ("geo_depends", Feel::po::value<std::string>()->default_value(""), "list of dependants file")
        ("geofile-path", Feel::po::value<std::string>()->default_value( "" ), "path to geofile ")
        ("repart", Feel::po::value<bool>()->default_value( false ), "repartition the mesh (use with care)")
        ("hsize", Feel::po::value<double>()->default_value( 0.5 ), "mesh size for coil")
        ("shape", Feel::po::value<std::string>()->default_value( "hypercube" ), "shape of the domain (either simplex or hypercube)")
        // Materials and Boundary conditions
        ("mat-bc-path", Feel::po::value<std::string>()->default_value( "." ), "path to access materials and boundary conditions")
        ("materials", Feel::po::value<std::string>(), "name of materials file")
        ("boundary-conditions", Feel::po::value<std::string>(), "name of boundary conditions file")
        // Mesh adaptation options
        ("meshadapt", Feel::po::value<bool>()->default_value( false ), "activate mesh adaptation")
        ("meshadapt_tol", Feel::po::value<double>()->default_value( 0.5 ), "tolerance parameter for mesh adaptation")
        ("meshadapt_method", Feel::po::value<std::string>()->default_value( "zz_map" ), "method for mesh adaptation: zz_map, hessian (L2 proj of solution) = hess1, hessian (L2 projection of hessian) = hess2")
        ("meshadapt_type", Feel::po::value<std::string>()->default_value( "isotropic" ), "type of mesh adaptation (isotropic, anisotropic)" )
        ("meshadapt_maxiter", Feel::po::value<int>()->default_value( 10 ), "max of mesh adaptation iterations")
        // Ginac so-file
        ("rebuild-so", Feel::po::value<bool>()->default_value( true ), "if true, remove all the so file (ginac) of the current directory")
        // Avoid export results
        ("export-results", Feel::po::value<bool>()->default_value( true ), "if false, results are not exported")
        ("export-highorders", Feel::po::value<bool>()->default_value( false ), "if false, high order results are not exported")
        //Units
        ("units", Feel::po::value<std::string>()->default_value( "m" ), "m (meters) or mm (millimeters)")
        // Common
        ("conductor_volume", Feel::po::value< std::vector<std::string> >()->default_value(std::vector<std::string>(), ""),
         "physical names of the volumes associated with conductor")
        ("conductor_js", Feel::po::value< std::vector<std::string> >()->default_value(std::vector<std::string>(), ""),
         "source current densities for associated conductor")
        ("insulator_volume", Feel::po::value< std::vector<std::string> >()->default_value(std::vector<std::string>(), ""),
         "physical names of the volumes associated with conductor")
        ("hchannels", Feel::po::value< std::vector<std::string> >()->default_value(std::vector<std::string>(), ""),
         "physical names of the boundary faces associated with longitudinal cooling channel")
        ("rchannels", Feel::po::value< std::vector<std::string> >()->default_value(std::vector<std::string>(), ""),
         "physical names of the boundary faces associated with radial cooling channel")        ;
    return hifimagnet_options
        .add( Feel::feel_options() );
}

