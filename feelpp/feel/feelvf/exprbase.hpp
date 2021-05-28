/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2014-03-01

  Copyright (C) 2014-2016 Feel++ Consortium

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
  \file exprbase.hpp
  \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
  \date 2010-04-21
*/
#ifndef FEELPP_EXPRBASE_HPP
#define FEELPP_EXPRBASE_HPP 1

#include <string>
#include <sstream>
#include <iostream>

#include <feel/feelpoly/context.hpp>

namespace Feel
{
/**
 * \class ExprBase
 * \brief Base class for expression
 *
 * Defines the common interface for all expression terms
 */
class ExprBase
{
public:


    /** @name Constants
     */
    //@{


    //@}

    /** @name Typedefs
     */
    //@{


    //@}

    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    ExprBase();

    //! destructor
    virtual ~ExprBase();

    //@}

    /** @name Operator overloads
     */
    //@{
    //@}

    /** @name Accessors
     */
    //@{

    /** Return a descriptive name for the expression subtype */
    virtual std::string typeName() const ;

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    /**
     * Write a simple text description suitable
     * for output to a terminal
     */
    virtual std::ostream& toText( std::ostream& os, bool paren ) const { return os; }

    /**
     * Write in a form suitable for LaTeX formatting
     */
    virtual std::ostream& toLatex( std::ostream& os, bool paren ) const { return os; }

    /**
     * write the expression into a std::string
     */
    std::string toString() const ;


    //@}



protected:

private:

};


namespace vf
{

class ExprDynamicBase
{
public:
    ExprDynamicBase() = default;
    explicit ExprDynamicBase( size_type c ) : M_context( c ) {}
    size_type dynamicContext() const { return M_context; }
private:
    size_type M_context = 0;
};

//!
//! @return true if the expr hass static context, false otherwise
//!
template <class T>
constexpr bool hasStaticContext()
{
    return !has_dynamic_v<T::context>;
}

//!
//! @return true if the expr hass static context, false otherwise
//!
template <class T>
inline bool hasStaticContext( T const& t )
{
    return !has_dynamic_v<T::context>;
}

//!
//! @return the static context  
//!
template <class T>
constexpr size_type staticContext()
{
    return T::context;
}

//!
//! @return the static context  
//!
template <class T>
inline size_type staticContext( T const& t )
{
    return T::context;
}

//!
//! @return true if the expression has dynamic context, false otherwise
//!
template <class T>
constexpr bool hasDynamicContext()
{
    return has_dynamic_v<T::context>;// && std::is_base_of_v<ExprDynamicBase, T>;
}

//!
//! @return true if the expression has dynamic context, false otherwise
//!
template <class T>
inline bool hasDynamicContext( T const& t )
{
    return has_dynamic_v<T::context>;// && std::is_base_of_v<ExprDynamicBase, T>;
}

//!
//! @return the dynamic context  if the expression has one or the static context otherwise
//!
template <class T>
size_type dynamicContext( T const& t )
{
    if constexpr ( hasDynamicContext<T>() )
    {
        return t.dynamicContext() | T::context;
    }
    else
    {
        return T::context;
    }
}

} // namespace vf

}
#endif /* FEELPP_EXPRBASE_HPP */
