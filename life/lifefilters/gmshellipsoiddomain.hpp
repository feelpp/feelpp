/* -*- mode: c++; coding: utf-8 -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-08-02

  Copyright (C) 2010 Université Joseph Fourier Grenoble 1

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
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file gmshellipsoiddomain.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2010-08-02
 */
#ifndef __GmshEllipsoidDomain_H
#define __GmshEllipsoidDomain_H 1


#include <boost/parameter.hpp>
#include <life/lifecore/life.hpp>
#include <life/lifecore/parameter.hpp>

#include <life/lifefilters/gmsh.hpp>
#include <life/lifemesh/simplex.hpp>

namespace Life
{
/**
 * \class GmshEllipsoidDomain
 *
 * \brief n-Ellipsoid Domain description for gmsh mesh generation.
 *
 *
 * \ingroup Importer
 * @author Christophe Prud'homme
 */
class GmshEllipsoidDomain : public Gmsh
{
    typedef Gmsh super;
public:


    /** @name Constants and Typedefs
     */
    //@{
    //@}

    /** @name Constructors, destructor
     */
    //@{

    GmshEllipsoidDomain( int Dim,
                         int Order= GMSH_ORDER_ONE,
                         DomainType dt = GMSH_REAL_DOMAIN );

    //! destructor
    ~GmshEllipsoidDomain() {}

    //@}


protected:

    std::string getDescription() const;

private:

    // 1D
    std::string getDescription1D() const;
    // 2D
    std::string getDescription2D() const;
    // 3D
    std::string getDescription3D() const;

};

} // Life


#endif /* __GmshEllipsoidDomain_H */
