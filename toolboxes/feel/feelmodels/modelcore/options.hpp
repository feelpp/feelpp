/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Vincent Chabannes <vincent.chabannes@feelpp.org.fr>
       Date: 2011-08-01

  Copyright (C) 2011 Universite Joseph Fourier (Grenoble I)

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
#ifndef FEELMODELS_OPTIONS_HPP
#define FEELMODELS_OPTIONS_HPP 1

#include <feel/options.hpp>

namespace Feel
{
Feel::po::options_description modelbase_options(std::string const& prefix);
Feel::po::options_description modelalgebraic_options(std::string const& prefix);
Feel::po::options_description modelnumerical_options(std::string const& prefix);
Feel::po::options_description densityviscosity_options(std::string const& prefix);
Feel::po::options_description fluidMechanics_options(std::string const& prefix="fluid");
Feel::po::options_description solidMechanics_options(std::string const& prefix="struct");
Feel::po::options_description alemesh_options(std::string const& prefix="fsi");
Feel::po::options_description fluidStructInteraction_options(std::string const& prefix);
Feel::po::options_description heat_options(std::string const& prefix="heat");
Feel::po::options_description electricity_options(std::string const& prefix);
Feel::po::options_description maxwell_options(std::string const& prefix);
Feel::po::options_description thermoElectric_options(std::string const& prefix);
Feel::po::options_description heatFluid_options(std::string const& prefix);
Feel::po::options_description advection_options(std::string const& prefix);
Feel::po::options_description redistanciation_fm_options(std::string const& prefix, bool addProjectorsOpts = true );
Feel::po::options_description redistanciation_hj_options(std::string const& prefix );
Feel::po::options_description levelset_options(std::string const& prefix);
Feel::po::options_description interfaceforces_options(std::string const& prefix);
Feel::po::options_description multifluid_options(std::string const& prefix, unsigned int nls = 3);
Feel::po::options_description coefficientformpde_options(std::string const& prefix);
Feel::po::options_description coefficientformpdes_options(std::string const& prefix);


Feel::po::options_description toolboxes_options(std::string const& type, std::string const& prefix);
inline
Feel::po::options_description toolboxes_options(std::string const& type) { return toolboxes_options( type,type ); }
}
#endif /* FEELMODELS_OPTIONS_HPP */
