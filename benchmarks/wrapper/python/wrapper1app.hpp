/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2010-06-09

  Copyright (C) 2010 Université Joseph Fourier (Grenoble I)

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
   \file wrapper1app.cpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2010-06-09
 */
//# marker1 #
#include <feel/options.hpp>
#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/application.hpp>
//# endmarker1 #

namespace Feel
{

/**
 * This routine defines some information about the application like
 * authors, version, or name of the application. The data returned is
 * typically used as an argument of a Feel::Application subclass.
 *
 * \return some data about the application.
 */
//# marker3 #
inline
AboutData
makeWrapper1About()
{
    AboutData about( "wrapper1" ,
                     "wrapper1" ,
                     "0.1",
                     "my Feel wrapper",
                     AboutData::License_GPL,
                     "Copyright (c) 2010 Universite Joseph Fourier" );

    about.addAuthor( "Christophe Prud'homme",
                     "developer",
                     "christophe.prudhomme@feelpp.org", "" );
    about.addAuthor( "Cyril Lamine",
                     "developer",
                     "", "" );
    return about;
}
//# endmarker3 #

/**
 * \class Wrapper1App
 *
 * This is a demo class to illustrate what is done (at the very least)
 * in subclasses of Feel::Application
 *
 */
//# marker4 #
class Wrapper1App: public Application
{
public:

    /**
     * constructor only about data and no options description
     */
    Wrapper1App( AboutData const& );

    void run( const double * X, unsigned long N,
              double * Y, unsigned long P )
    {

        Y[0] = X[0];
    }

};
//# endmarker4 #

//# marker5 #
Wrapper1App::Wrapper1App( AboutData const& ad )
    :
    Application( ad )
{}
//# endmarker5 #
}
