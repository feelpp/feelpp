/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2005-07-17

  Copyright (C) 2005,2006 EPFL
  Copyright (C) 2007 Universite Joseph Fourier (Grenoble I)

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
   \file binaryfunctor.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-07-17
 */
#ifndef __BinaryFunctor_H
#define __BinaryFunctor_H 1

#include <feel/feelvf/functordomain.hpp>

namespace Feel
{
namespace vf
{
/// \cond detail
/**
 * \class BinaryFunctor
 * \brief unary functors base class
 *
 * @author Christophe Prud'homme
 * @see
 */
template<typename T1 = double,typename T2 = double>
class BinaryFunctor
{
public:

    typedef T1 value_1_type;
    typedef T2 value_2_type;
    typedef FunctorDomain<value_1_type> functordomain_1_type;
    typedef FunctorDomain<value_2_type> functordomain_2_type;
    //typedef boost::shared_ptr<functordomain_type> functordomain_ptrtype;
    /** @name Constructors, destructor
     */
    //@{

    BinaryFunctor( std::string const& name,
                   std::pair<functordomain_1_type,functordomain_2_type> const& domain = std::make_pair( UnboundedDomain<value_1_type>(),
                                                                                                        UnboundedDomain<value_2_type>() ) )
        :
        M_name( name ),
        M_domain( domain )
    {}

    BinaryFunctor( BinaryFunctor const & uf )
        :
        M_name( uf.M_name ),
        M_domain( uf.M_domain )
    {}

    virtual ~BinaryFunctor()
    {}

    //@}

    /** @name Accessors
     */
    //@{

    /**
     * @return functor name
     */
    std::string const& name() const
    {
        return M_name;
    }


    /**
     * @return functor domain object
     */
    std::pair<functordomain_1_type,functordomain_2_type> const& domain() const
    {
        return M_domain;
    }

    //@}

    /** @name  Methods
     */
    //@{

    /**
     * evaluate the functor at \c nx points \c x and store the values
     *  in the array \c f
     *
     * \param nx number of points for function evaluation
     * \param x  points for function evaluation
     * \param f  store function evaluation at nx points
     */
    virtual void eval( int nx, value_1_type const* x, value_2_type const* y, value_1_type* f ) const = 0;


    //@}



protected:

private:
    std::string M_name;

    std::pair<functordomain_1_type,functordomain_2_type> M_domain;

};
/// \endcond
}
}
#endif /* __BinaryFunctor_H */
