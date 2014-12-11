/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-07-21

  Copyright (C) 2007 Universite Joseph Fourier (Grenoble I)

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
/**
   \file options.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-07-21
 */
#ifndef FEELPP_OPTIONS_HPP
#define FEELPP_OPTIONS_HPP 1

#include <boost/program_options.hpp>

namespace Feel
{
namespace po = boost::program_options;


po::options_description
file_options( std::string const& prefix );

po::options_description
generic_options();

po::options_description
feel_options( std::string const& prefix = "" );

inline po::options_description
feel_nooptions() { return po::options_description(); }

//! add benchmark options to feel++ applications
po::options_description benchmark_options( std::string const& prefix = "" );

/**
 * \param prefix prefix given to the  backend option
 * \return backend command line options description
 */
po::options_description backend_options( std::string const& prefix = "" );

po::options_description backendpetsc_options( std::string const& prefix = "" );

po::options_description mesh_options( int Dim, std::string const& prefix = "" );

po::options_description gmsh_options( std::string const& prefix = "" );

po::options_description gmsh_domain_options( std::string const& prefix = "" );

/**
 * command line options for Onelab interface
 */
po::options_description onelab_options( std::string const& prefix = "" );

/**
 * command line options for multithreading and GPU parallelization
 */
po::options_description parallel_options( std::string const& prefix = "" );

po::options_description ginac_options( std::string const& prefix = "" );

/**
 * defines solver eigen options
 *
 * The \p prefix parameter allows to set different eigensolver options for
 * different eigensolver. It allows to distinguish between these options
 * \code
 * // register two slepc eigensolver options
 * add_options( solvereigen_options( "eigen1" ) ).add_options( solvereigen_options( "eigen2" ));
 * // build an eigen solver associated with option set eigen1
 * SolverEigen<double>::build( vm, "eigen1" );
 * \endcode
 *
 * \param prefix prefix allows to prefix options
 */
po::options_description solvereigen_options( std::string const& prefix = "" );

/**
 * command line options for non linear solver
 */
po::options_description nlsolver_options();

/**
 * command line options for BDF
 */
po::options_description on_options( std::string const& prefix = "" );

po::options_description ts_options( std::string const& prefix = "" );
po::options_description bdf_options( std::string const& prefix = "" );

/**
 * command line options for exporter
 */

po::options_description exporter_options( std::string const& prefix = "" );

po::options_description material_options( std::string const& prefix = "" );

po::options_description error_options( std::string const& prefix = "" );

po::options_description functionspace_options( std::string const& prefix = "" );

po::options_description aitken_options( std::string const& prefix = "" );

}
#endif // FEELPP_OPTIONS_HPP
