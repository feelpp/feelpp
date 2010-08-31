/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
   \author Christophe Prud'homme <christophe.prudhomme@ujf-grenoble.fr>
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
template<int Dim, int Order, int RDim, template<uint16_type, uint16_type, uint16_type> class Entity >
class GmshHypercubeDomain : public Gmsh
{
    typedef Gmsh super;
public:


    /** @name Constants and Typedefs
     */
    //@{
    static const uint16_type nDim = Dim;
    static const uint16_type nOrder = Order;
    static const uint16_type nRealDim = RDim;

    typedef Entity<Dim,Order, nRealDim> entity_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    GmshHypercubeDomain()
        :
        super(Dim, Order)
    {
    }

    GmshHypercubeDomain( GmshHypercubeDomain const & td )
        :
        super( td )
    {
    }
    ~GmshHypercubeDomain()
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
    {return getDescription( mpl::int_<nDim>(), mpl::bool_<entity_type::is_simplex_product>() );}
    // 1D
    std::string getDescription( mpl::int_<1>, mpl::bool_<false> ) const;
    std::string getDescription( mpl::int_<1>, mpl::bool_<true> ) const
    { return getDescription( mpl::int_<1>(), mpl::bool_<false>() ); }
    // 2D
    std::string getDescription( mpl::int_<2>, mpl::bool_<false> ) const;
    std::string getDescription( mpl::int_<2>, mpl::bool_<true> ) const;
    // 3D
    std::string getDescription( mpl::int_<3>, mpl::bool_<false>, bool do_recombine = false ) const;
    std::string getDescription( mpl::int_<3>, mpl::bool_<true> ) const;

private:


};

} // Feel

#if !defined( FEEL_INSTANTIATION_MODE )
# include <feel/feelfilters/gmshhypercubedomain.cpp>
#endif // FEEL_INSTANTIATION_MODE

#endif /* __GmshHypercubeDomain_H */
