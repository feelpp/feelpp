/* -*- mode: c++; coding: utf-8 -*-

  This file is part of the Life library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
       Date: 2007-10-11

  Copyright (C) 2007 Université Joseph Fourier Grenoble 1

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
   \file gmshsimplexdomain.hpp
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
   \date 2007-10-11
 */
#ifndef __GmshSimplexDomain_H
#define __GmshSimplexDomain_H 1


#include <boost/parameter.hpp>
#include <life/lifecore/life.hpp>
#include <life/lifecore/parameter.hpp>

#include <life/lifefilters/gmsh.hpp>
#include <life/lifemesh/simplex.hpp>

namespace Life
{
/**
 * \class GmshSimplexDomain
 * \brief Simplex Domain description for gmsh mesh generation
 *
 * \ingroup Importer
 * @author Christophe Prud'homme
 */
template<int Dim, int Order>
class GmshSimplexDomain : public Gmsh
{
    typedef Gmsh super;
public:


    /** @name Constants and Typedefs
     */
    //@{

    static const uint16_type nDim = Dim;
    static const uint16_type nOrder = Order;

    typedef Simplex<Dim,Order, Dim> entity_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    GmshSimplexDomain( DomainType dt = GMSH_REAL_DOMAIN );

    GmshSimplexDomain( GmshSimplexDomain const & td )
        :
        super( td ),
        _M_descr( td._M_descr )
    {
    }
    ~GmshSimplexDomain()
    {}

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
    std::string getDescription() const
    {return getDescription( mpl::int_<nDim>() );}
    // 1D
    std::string getDescription( mpl::int_<1> ) const;
    // 2D
    std::string getDescription( mpl::int_<2> ) const;
    // 3D
    std::string getDescription( mpl::int_<3> ) const;

private:

    std::string _M_descr;
};

} // Life

#if !defined( LIFE_INSTANTIATION_MODE )
# include <life/lifefilters/gmshsimplexdomain.cpp>
#endif // LIFE_INSTANTIATION_MODE

#endif /* __GmshSimplexDomain_H */
