/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2006-11-26

  Copyright (C) 2006 Université Joseph Fourier

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
   \file gmshhypercubedomain.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2006-11-26
 */
#ifndef __GmshHypercubeDomain_H
#define __GmshHypercubeDomain_H 1

#include <feel/feelfilters/gmsh.hpp>

namespace Feel
{
/**
 * \class GmshHypercubeDomain
 * \brief Tensorized Domain description for gmsh mesh generation
 *
 * \ingroup Importer
 * @author Christophe Prud'homme
 */
class GmshHypercubeDomain : public Gmsh
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

    GmshHypercubeDomain( int dim, int order );
    GmshHypercubeDomain( int dim, int order, int rdim, bool use_hypercube );
    GmshHypercubeDomain( GmshHypercubeDomain const & td );
    ~GmshHypercubeDomain();

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

    //@}



private:
    int M_rdim;
    bool M_use_hypercube;
    //
    std::string getDescription() const;
    // 1D
    std::string getDescription1D() const;
    // 2D
    std::string getDescription2D() const;
    // 3D
    std::string getDescription3D() const;

private:


};

} // Feel

#endif /* __GmshHypercubeDomain_H */
