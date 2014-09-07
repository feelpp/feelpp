/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2010-04-27

  Copyright (C) 2010-2013 Feel++ Consortium

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
   \file formcontextbase.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2010-04-27
 */
#ifndef FEELPP_FORMCONTEXTBASE_HPP
#define FEELPP_FORMCONTEXTBASE_HPP 1

#include <feel/feelvf/detail/gmc.hpp>

namespace Feel
{
namespace vf
{
/// \cond detail
namespace detail
{

/**
 * \class FormContextBase
 * \brief base class for bi/linear form contexts
 */
template<typename GeomapContext, typename IM, typename GeomapExprContext>
class FormContextBase
{
public:

    /** @name Typedefs
     */
    //@{
    typedef GeomapContext map_geometric_mapping_context_type;
    typedef typename fusion::result_of::value_at_key<GeomapContext,gmc<0> >::type geometric_mapping_context_ptrtype;
    typedef typename geometric_mapping_context_ptrtype::element_type geometric_mapping_context_type;
    typedef typename geometric_mapping_context_type::gm_type geometric_mapping_type;

    typedef mpl::int_<fusion::result_of::template size<GeomapContext>::type::value> map_size;

    typedef typename mpl::if_<mpl::equal_to<map_size,mpl::int_<2> >, gmc<1>, gmc<0> >::type gmc1;

    typedef typename fusion::result_of::value_at_key<GeomapContext,gmc<0> >::type left_gmc_ptrtype;
    typedef typename fusion::result_of::value_at_key<GeomapContext,gmc<0> >::type::element_type left_gmc_type;
    typedef typename fusion::result_of::value_at_key<GeomapContext,gmc1 >::type right_gmc_ptrtype;
    typedef typename fusion::result_of::value_at_key<GeomapContext,gmc1 >::type::element_type right_gmc_type;

    typedef fusion::map<fusion::pair<gmc<0>, left_gmc_ptrtype> > map_left_gmc_type;
    typedef fusion::map<fusion::pair<gmc<0>, right_gmc_ptrtype> > map_right_gmc_type;


    typedef GeomapExprContext map_geometric_mapping_expr_context_type;

    //@}

    virtual ~FormContextBase() {}

    virtual void update( map_geometric_mapping_context_type const& _gmc,
                         map_geometric_mapping_expr_context_type const& _gmcExpr ) = 0;

    virtual void update( map_geometric_mapping_context_type const& _gmc,
                         map_geometric_mapping_expr_context_type const& _gmcExpr,
                         mpl::int_<2> ) = 0;


    virtual void update( map_geometric_mapping_context_type const& _gmc,
                         map_geometric_mapping_expr_context_type const& _gmcExpr,
                         IM const& im ) = 0;

    virtual void update( map_geometric_mapping_context_type const& _gmc,
                         map_geometric_mapping_expr_context_type const& _gmcExpr,
                         IM const& im, mpl::int_<2> ) = 0;


    virtual void integrate() = 0;

    virtual void assemble() = 0;

    virtual void assemble( mpl::int_<2> ) = 0;
};
}
/// \endcond
}
}
#endif /* __FormContextBase_H */
