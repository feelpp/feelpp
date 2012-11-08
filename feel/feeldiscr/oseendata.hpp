/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-11-20

  Copyright (C) 2008-2011 Université Joseph Fourier (Grenoble I)

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
   \file oseendata.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-11-20
 */
#ifndef __OseenData_H
#define __OseenData_H 1

namespace Feel
{
/*!
  \class OseenDefaults
  \brief Oseen Default Data

  @author Christophe Prud'homme
  @see
*/
// struct to hold default values for Oseen
struct OseenDefaults
{
    // Default constructor with default default values
    OseenDefaults()
        :
        BC_COEFF_DIFF( 100.0 ),
        BC_COEFF_CONV( 100.0 ),
        STAB_COEFF_DIV( 0.0 ),
        STAB_COEFF_P( 0.0 ),
        EPS_COMPRESS( 0.0 ),
        DIVDIV_COEFF( 0.0 ),
        WEAK_DIRICHLET( true ),
        EXPORT_MATLAB( false )
    {}

    // coefficient for diffusive terms of weak Dirichlet conditions
    double BC_COEFF_DIFF;
    // coefficient for convective terms of weak Dirichlet conditions
    double BC_COEFF_CONV;
    // coefficient for divergence jump stabilisation
    double STAB_COEFF_DIV;
    // coefficient for pressure gradient jump stabilization
    double STAB_COEFF_P;
    // coefficient for pseudo-compressibility term
    double EPS_COMPRESS;
    // coefficient for divergence penalty term
    double DIVDIV_COEFF;
    // whether to use weak instead of strong dirichlet boundary conditions
    bool WEAK_DIRICHLET;
    // whether to export matrix and vector in matlab format
    bool EXPORT_MATLAB;
};

po::options_description oseen_options( std::string const& prefix="",
                                       OseenDefaults defaults=OseenDefaults() );

} // Feel
#endif /* __OseenData_H */
