/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
       Date: 2006-11-16

  Copyright (C) 2006 EPFL

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
   \file operator.hpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2006-11-16
 */

#ifndef _OPERATOR_HPP_
#define _OPERATOR_HPP_

#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/fsfunctionallinear.hpp>

namespace Feel
{

/**
 * \class Operator
 * \brief Operator between function spaces
 */
template<class DomainSpace, class DualImageSpace>
class Operator
{
public:

    // -- TYPEDEFS --
    typedef DomainSpace     domain_space_type;
    typedef DualImageSpace  dual_image_space_type;

    typedef typename domain_space_type::value_type value_type;
    typedef Backend<value_type> backend_type;
    typedef boost::shared_ptr<backend_type> backend_ptrtype;

    typedef boost::shared_ptr<domain_space_type> domain_space_ptrtype;
    typedef boost::shared_ptr<dual_image_space_type>
    dual_image_space_ptrtype;
    typedef typename domain_space_type::element_type domain_element_type;
    typedef FsFunctionalLinear<dual_image_space_type> image_element_type;

    Operator()
        :
        M_domainSpace(),
        M_dualImageSpace()
    {}

    Operator( domain_space_ptrtype     domainSpace,
              dual_image_space_ptrtype dualImageSpace ) :
        M_domainSpace( domainSpace ),
        M_dualImageSpace( dualImageSpace )
    {}

    virtual ~Operator() {}

    void setDomainSpace( domain_space_ptrtype const& domainspace )
    {
        M_domainSpace = domainspace;
    }

    void setDualImageSpace( dual_image_space_ptrtype const& dualImageSpace )
    {
        M_dualImageSpace = dualImageSpace;
    }

    // apply the operator: ie := Op de
    virtual void
    apply( const domain_element_type& de,
           image_element_type&        ie ) const = 0;



    // for convenience
    image_element_type operator()( const domain_element_type& de ) const
    {
        image_element_type ie( M_dualImageSpace );
        apply( de, ie );
        return ie;
    }

    domain_space_ptrtype domainSpace()
    {
        return M_domainSpace;
    }

    dual_image_space_ptrtype dualImageSpace()
    {
        return M_dualImageSpace;
    }

    const domain_space_ptrtype domainSpace() const
    {
        return M_domainSpace;
    }

    const dual_image_space_ptrtype dualImageSpace() const
    {
        return M_dualImageSpace;
    }

private:

    domain_space_ptrtype     M_domainSpace;
    dual_image_space_ptrtype M_dualImageSpace;

}; // class Operator

} // Feel

#endif /* _OPERATOR_HPP_ */
