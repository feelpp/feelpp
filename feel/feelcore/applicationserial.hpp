/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-10-18

  Copyright (C) 2005,2006,2009 EPFL

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
   \file application.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-10-18
 */
#ifndef __Application_H
#define __Application_H 1

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <feel/feelconfig.h>

#include <feel/feelcore/feel.hpp>
#include <feel/feelcore/application.hpp>

namespace Feel
{
/**
 * \class Application
 *\ingroup Core
 *\brief SERIAL Application
 *
 * @author Christophe Prud'homme
 * @see
 */
class Application : public Application
{
    typedef Application super;
public:


    /** @name Typedefs
     */
    //@{


    //@}

    /** @name Constructors, destructor
     */
    //@{

    /**
     * Construct an SERIAL Application
     *
     * @param argc number of arguments on the command line
     * @param argv arguments in the command line
     * @param ad \p AboutData structure for this \p Application
     */
    Application( int argc, char** argv, AboutData const& ad );

    /**
     * Construct an SERIAL Application
     *
     * @param argc number of arguments on the command line
     * @param argv arguments in the command line
     * @param ad \p AboutData structure for this \p Application
     * @param od \p po::options_description structure for this \p Application
     */
    Application( int argc, char** argv, AboutData const& ad, po::options_description const& od );

    Application( Application const & );

    ~Application();

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    template<class T>
    static void Broadcast( T& /*obj*/, int /*root*/ = 0 )
    {}

    /**
     * @return the barrier
     */
    static void barrier() { }

    //@}



protected:


private:



private:

};


}
#endif /* __Application_H */
