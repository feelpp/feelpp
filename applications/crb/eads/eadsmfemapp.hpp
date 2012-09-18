/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2010-02-06

  Copyright (C) 2010-2011 Université Joseph Fourier (Grenoble I)

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
 *  @file  eadsmfem.hpp
 *  @brief The header file for opuseadsmfem
 */
#ifndef OPUSEADSMFEM_H
#define OPUSEADSMFEM_H

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/environment.hpp>
#include <feel/feelcore/application.hpp>
#include <feel/options.hpp>

#include <eads.hpp>

#include <feel/feelcrb/crbwrapper.hpp>
#include <opusdata.hpp>
#include <opusmodelbase.hpp>
#include <opusmodelfactory.hpp>


/*
 *  This is the declaration of function named 'opuseadsmfem' into the wrapper.
 */
namespace Feel
{
/**
 * \class EadsMFemApp
 * \brief EadsMFem application
 *
 * This class implements the Eads MFem application, getting the command line
 * arguments and running the actual code.
 *
 * @author Christophe Prud'homme
 */
class EadsMFemApp   : public Application
{
    typedef Application super;
public:

    typedef OpusModelBase opus_type;
    typedef boost::shared_ptr<opus_type> opus_ptrtype;

    EadsMFemApp( AboutData const& ad, po::options_description const& od )
        :
        super( ad, od )
    {
        this->init();
    }
    EadsMFemApp( int argc, char** argv, AboutData const& ad, po::options_description const& od )
        :
        super( argc, argv, ad, od )
    {
        this->init();
    }
    void init();
    void run();
    void run( const double * X, unsigned long N,
              double * Y, unsigned long P );

    opus_ptrtype opus;

}; // Opus
} // Feel


#endif /* OPUSEADSMFEM_H */

