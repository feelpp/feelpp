/*
   adfunctions.hpp	ADType expression templates
   This file is part of gstlibs.

   Copyright (C) 2011 Christophe Prud'homme

   gstlibs is free software; you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   gstlibs is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with gstlibs; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

// Generated source file.  Do not edit.
// /home/prudhomm/Devel/FEEL/feel/feel/feelopt/adgenerate.cpp Jun 22 2011 01:46:56

#ifndef AD_FUNCS_HPP
#define AD_FUNCS_HPP

namespace Feel
{
/****************************************************************************
 * Cosinus function
 ****************************************************************************/
/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-02-14

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

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
   \file adfunctions.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-14
 */
template <class Expr> class ADFuncCos
{
public:
    enum { nvar = Expr::nvar };
    typedef typename Expr::value_type value_type;
protected:
    ADFuncCos () {}

    Expr expr_;
public:

    ADFuncCos ( const Expr & expr ) : expr_( expr )
    {
        ;
    }


    template<int isFundamental, typename Expr_>
    struct Value
    {

        typedef typename Expr_::value_type value_type;
        static value_type signum( value_type a )
        {
            if ( abs( a ) < std::numeric_limits<value_type>::epsilon() )
                return value_type( 0 );

            return a/abs( a );
        }
        static value_type value( Expr_ const& expr_ )
        {
            return cos( expr_.value() );
        }
        static value_type grad( Expr_ const& expr_, int __i )
        {
            return -expr_.grad( __i )*sin( expr_.value() );
        }
        static value_type hessian( Expr_ const& expr_, int __i, int __j )
        {
            return -expr_.hessian( __i,__j )*sin( expr_.value() ) - expr_.grad( __i )*expr_.grad( __j ) * cos( expr_.value() );
        }
    };
    template<typename Expr_>
    struct Value<true, Expr_>
    {
        typedef typename Expr_::value_type value_type;
        static value_type signum( value_type a )
        {
            if ( std::abs( a ) < std::numeric_limits<value_type>::epsilon() )
                return value_type( 0 );

            return a/std::abs( a );
        }
        static value_type value( Expr_ const& expr_ )
        {
            return std::cos( expr_.value() );
        }
        static value_type grad( Expr_ const& expr_, int __i )
        {
            return -expr_.grad( __i )*std::sin( expr_.value() );
        }
        static value_type hessian( Expr_ const& expr_, int __i, int __j )
        {
            return -expr_.hessian( __i,__j )*std::sin( expr_.value() ) - expr_.grad( __i )*expr_.grad( __j ) * std::cos( expr_.value() );
        }
    };

    inline value_type value() const
    {
        return Value<true/*boost::type_traits::is_fundamental<value_type>::value*/,Expr>::value( expr_ );
    }
    inline value_type grad( int __i ) const
    {
        return  Value<true/*boost::type_traits::is_fundamental<value_type>::value*/,Expr>::grad( expr_, __i );
    }
    inline value_type hessian( int __i, int __j ) const
    {
        return  Value<true/*boost::type_traits::is_fundamental<value_type>::value*/,Expr>::hessian( expr_, __i, __j );
    }
    inline bool deps( int __i ) const
    {
        return expr_.deps( __i ) ;
    }
};

template <class Expr> inline ADExpr< ADFuncCos< ADExpr<Expr> > >
cos ( const ADExpr<Expr>& expr )
{
    typedef ADFuncCos< ADExpr<Expr> > expr_t;
    return ADExpr< expr_t >(  expr_t( expr ) );
}

template <class T, int Nvar, int Order, int Var> inline ADExpr< ADFuncCos< ADType<T, Nvar, Order, Var> > >
cos ( const ADType<T, Nvar, Order, Var>& x )
{
    typedef ADFuncCos< ADType<T, Nvar, Order, Var> > expr_t;
    return ADExpr< expr_t >(  expr_t( x ) );
}
/****************************************************************************
 * Sinus function
 ****************************************************************************/
/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-02-14

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

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
   \file adfunctions.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-14
 */
template <class Expr> class ADFuncSin
{
public:
    enum { nvar = Expr::nvar };
    typedef typename Expr::value_type value_type;
protected:
    ADFuncSin () {}

    Expr expr_;
public:

    ADFuncSin ( const Expr & expr ) : expr_( expr )
    {
        ;
    }


    template<int isFundamental, typename Expr_>
    struct Value
    {

        typedef typename Expr_::value_type value_type;
        static value_type signum( value_type a )
        {
            if ( abs( a ) < std::numeric_limits<value_type>::epsilon() )
                return value_type( 0 );

            return a/abs( a );
        }
        static value_type value( Expr_ const& expr_ )
        {
            return sin( expr_.value() );
        }
        static value_type grad( Expr_ const& expr_, int __i )
        {
            return expr_.grad( __i )*cos( expr_.value() );
        }
        static value_type hessian( Expr_ const& expr_, int __i, int __j )
        {
            return expr_.hessian( __i,__j )*cos( expr_.value() ) - expr_.grad( __i )*expr_.grad( __j ) * sin( expr_.value() );
        }
    };
    template<typename Expr_>
    struct Value<true, Expr_>
    {
        typedef typename Expr_::value_type value_type;
        static value_type signum( value_type a )
        {
            if ( std::abs( a ) < std::numeric_limits<value_type>::epsilon() )
                return value_type( 0 );

            return a/std::abs( a );
        }
        static value_type value( Expr_ const& expr_ )
        {
            return std::sin( expr_.value() );
        }
        static value_type grad( Expr_ const& expr_, int __i )
        {
            return expr_.grad( __i )*std::cos( expr_.value() );
        }
        static value_type hessian( Expr_ const& expr_, int __i, int __j )
        {
            return expr_.hessian( __i,__j )*std::cos( expr_.value() ) - expr_.grad( __i )*expr_.grad( __j ) * std::sin( expr_.value() );
        }
    };

    inline value_type value() const
    {
        return Value<true/*boost::type_traits::is_fundamental<value_type>::value*/,Expr>::value( expr_ );
    }
    inline value_type grad( int __i ) const
    {
        return  Value<true/*boost::type_traits::is_fundamental<value_type>::value*/,Expr>::grad( expr_, __i );
    }
    inline value_type hessian( int __i, int __j ) const
    {
        return  Value<true/*boost::type_traits::is_fundamental<value_type>::value*/,Expr>::hessian( expr_, __i, __j );
    }
    inline bool deps( int __i ) const
    {
        return expr_.deps( __i ) ;
    }
};

template <class Expr> inline ADExpr< ADFuncSin< ADExpr<Expr> > >
sin ( const ADExpr<Expr>& expr )
{
    typedef ADFuncSin< ADExpr<Expr> > expr_t;
    return ADExpr< expr_t >(  expr_t( expr ) );
}

template <class T, int Nvar, int Order, int Var> inline ADExpr< ADFuncSin< ADType<T, Nvar, Order, Var> > >
sin ( const ADType<T, Nvar, Order, Var>& x )
{
    typedef ADFuncSin< ADType<T, Nvar, Order, Var> > expr_t;
    return ADExpr< expr_t >(  expr_t( x ) );
}
/****************************************************************************
 * Tangent function
 ****************************************************************************/
/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-02-14

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

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
   \file adfunctions.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-14
 */
template <class Expr> class ADFuncTan
{
public:
    enum { nvar = Expr::nvar };
    typedef typename Expr::value_type value_type;
protected:
    ADFuncTan () {}

    Expr expr_;
public:

    ADFuncTan ( const Expr & expr ) : expr_( expr )
    {
        ;
    }


    template<int isFundamental, typename Expr_>
    struct Value
    {

        typedef typename Expr_::value_type value_type;
        static value_type signum( value_type a )
        {
            if ( abs( a ) < std::numeric_limits<value_type>::epsilon() )
                return value_type( 0 );

            return a/abs( a );
        }
        static value_type value( Expr_ const& expr_ )
        {
            return tan( expr_.value() );
        }
        static value_type grad( Expr_ const& expr_, int __i )
        {
            return expr_.grad( __i )*( 1.+tan( expr_.value() )*tan( expr_.value() ) );
        }
        static value_type hessian( Expr_ const& expr_, int __i, int __j )
        {
            return 2*tan( expr_.value() )*( 1+tan( expr_.value() )*tan( expr_.value() ) )*expr_.grad( __i )*expr_.grad( __j ) + ( 1+tan( expr_.value() )*tan( expr_.value() ) )*expr_.hessian( __i,__j );
        }
    };
    template<typename Expr_>
    struct Value<true, Expr_>
    {
        typedef typename Expr_::value_type value_type;
        static value_type signum( value_type a )
        {
            if ( std::abs( a ) < std::numeric_limits<value_type>::epsilon() )
                return value_type( 0 );

            return a/std::abs( a );
        }
        static value_type value( Expr_ const& expr_ )
        {
            return std::tan( expr_.value() );
        }
        static value_type grad( Expr_ const& expr_, int __i )
        {
            return expr_.grad( __i )*( 1.+std::tan( expr_.value() )*std::tan( expr_.value() ) );
        }
        static value_type hessian( Expr_ const& expr_, int __i, int __j )
        {
            return 2*std::tan( expr_.value() )*( 1+std::tan( expr_.value() )*std::tan( expr_.value() ) )*expr_.grad( __i )*expr_.grad( __j ) + ( 1+std::tan( expr_.value() )*std::tan( expr_.value() ) )*expr_.hessian( __i,__j );
        }
    };

    inline value_type value() const
    {
        return Value<true/*boost::type_traits::is_fundamental<value_type>::value*/,Expr>::value( expr_ );
    }
    inline value_type grad( int __i ) const
    {
        return  Value<true/*boost::type_traits::is_fundamental<value_type>::value*/,Expr>::grad( expr_, __i );
    }
    inline value_type hessian( int __i, int __j ) const
    {
        return  Value<true/*boost::type_traits::is_fundamental<value_type>::value*/,Expr>::hessian( expr_, __i, __j );
    }
    inline bool deps( int __i ) const
    {
        return expr_.deps( __i ) ;
    }
};

template <class Expr> inline ADExpr< ADFuncTan< ADExpr<Expr> > >
tan ( const ADExpr<Expr>& expr )
{
    typedef ADFuncTan< ADExpr<Expr> > expr_t;
    return ADExpr< expr_t >(  expr_t( expr ) );
}

template <class T, int Nvar, int Order, int Var> inline ADExpr< ADFuncTan< ADType<T, Nvar, Order, Var> > >
tan ( const ADType<T, Nvar, Order, Var>& x )
{
    typedef ADFuncTan< ADType<T, Nvar, Order, Var> > expr_t;
    return ADExpr< expr_t >(  expr_t( x ) );
}
/****************************************************************************
 * Sqrt function
 ****************************************************************************/
/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-02-14

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

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
   \file adfunctions.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-14
 */
template <class Expr> class ADFuncSqrt
{
public:
    enum { nvar = Expr::nvar };
    typedef typename Expr::value_type value_type;
protected:
    ADFuncSqrt () {}

    Expr expr_;
public:

    ADFuncSqrt ( const Expr & expr ) : expr_( expr )
    {
        ;
    }


    template<int isFundamental, typename Expr_>
    struct Value
    {

        typedef typename Expr_::value_type value_type;
        static value_type signum( value_type a )
        {
            if ( abs( a ) < std::numeric_limits<value_type>::epsilon() )
                return value_type( 0 );

            return a/abs( a );
        }
        static value_type value( Expr_ const& expr_ )
        {
            return sqrt( expr_.value() );
        }
        static value_type grad( Expr_ const& expr_, int __i )
        {
            return expr_.grad( __i )/( 2.*sqrt( expr_.value() ) );
        }
        static value_type hessian( Expr_ const& expr_, int __i, int __j )
        {
            return expr_.hessian( __i,__j )/( 2.*sqrt( expr_.value() ) )-expr_.grad( __i )*expr_.grad( __j )/( 4.*pow( expr_.value(),1.5 ) );
        }
    };
    template<typename Expr_>
    struct Value<true, Expr_>
    {
        typedef typename Expr_::value_type value_type;
        static value_type signum( value_type a )
        {
            if ( std::abs( a ) < std::numeric_limits<value_type>::epsilon() )
                return value_type( 0 );

            return a/std::abs( a );
        }
        static value_type value( Expr_ const& expr_ )
        {
            return std::sqrt( expr_.value() );
        }
        static value_type grad( Expr_ const& expr_, int __i )
        {
            return expr_.grad( __i )/( 2.*std::sqrt( expr_.value() ) );
        }
        static value_type hessian( Expr_ const& expr_, int __i, int __j )
        {
            return expr_.hessian( __i,__j )/( 2.*std::sqrt( expr_.value() ) )-expr_.grad( __i )*expr_.grad( __j )/( 4.*std::pow( expr_.value(),1.5 ) );
        }
    };

    inline value_type value() const
    {
        return Value<true/*boost::type_traits::is_fundamental<value_type>::value*/,Expr>::value( expr_ );
    }
    inline value_type grad( int __i ) const
    {
        return  Value<true/*boost::type_traits::is_fundamental<value_type>::value*/,Expr>::grad( expr_, __i );
    }
    inline value_type hessian( int __i, int __j ) const
    {
        return  Value<true/*boost::type_traits::is_fundamental<value_type>::value*/,Expr>::hessian( expr_, __i, __j );
    }
    inline bool deps( int __i ) const
    {
        return expr_.deps( __i ) ;
    }
};

template <class Expr> inline ADExpr< ADFuncSqrt< ADExpr<Expr> > >
sqrt ( const ADExpr<Expr>& expr )
{
    typedef ADFuncSqrt< ADExpr<Expr> > expr_t;
    return ADExpr< expr_t >(  expr_t( expr ) );
}

template <class T, int Nvar, int Order, int Var> inline ADExpr< ADFuncSqrt< ADType<T, Nvar, Order, Var> > >
sqrt ( const ADType<T, Nvar, Order, Var>& x )
{
    typedef ADFuncSqrt< ADType<T, Nvar, Order, Var> > expr_t;
    return ADExpr< expr_t >(  expr_t( x ) );
}
/****************************************************************************
 * Exponential function
 ****************************************************************************/
/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-02-14

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

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
   \file adfunctions.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-14
 */
template <class Expr> class AdFuncExp
{
public:
    enum { nvar = Expr::nvar };
    typedef typename Expr::value_type value_type;
protected:
    AdFuncExp () {}

    Expr expr_;
public:

    AdFuncExp ( const Expr & expr ) : expr_( expr )
    {
        ;
    }


    template<int isFundamental, typename Expr_>
    struct Value
    {

        typedef typename Expr_::value_type value_type;
        static value_type signum( value_type a )
        {
            if ( abs( a ) < std::numeric_limits<value_type>::epsilon() )
                return value_type( 0 );

            return a/abs( a );
        }
        static value_type value( Expr_ const& expr_ )
        {
            return exp( expr_.value() );
        }
        static value_type grad( Expr_ const& expr_, int __i )
        {
            return expr_.grad( __i )*exp( expr_.value() );
        }
        static value_type hessian( Expr_ const& expr_, int __i, int __j )
        {
            return expr_.hessian( __i,__j )*exp( expr_.value() ) + expr_.grad( __i )*expr_.grad( __j )*exp( expr_.value() );
        }
    };
    template<typename Expr_>
    struct Value<true, Expr_>
    {
        typedef typename Expr_::value_type value_type;
        static value_type signum( value_type a )
        {
            if ( std::abs( a ) < std::numeric_limits<value_type>::epsilon() )
                return value_type( 0 );

            return a/std::abs( a );
        }
        static value_type value( Expr_ const& expr_ )
        {
            return std::exp( expr_.value() );
        }
        static value_type grad( Expr_ const& expr_, int __i )
        {
            return expr_.grad( __i )*std::exp( expr_.value() );
        }
        static value_type hessian( Expr_ const& expr_, int __i, int __j )
        {
            return expr_.hessian( __i,__j )*std::exp( expr_.value() ) + expr_.grad( __i )*expr_.grad( __j )*std::exp( expr_.value() );
        }
    };

    inline value_type value() const
    {
        return Value<true/*boost::type_traits::is_fundamental<value_type>::value*/,Expr>::value( expr_ );
    }
    inline value_type grad( int __i ) const
    {
        return  Value<true/*boost::type_traits::is_fundamental<value_type>::value*/,Expr>::grad( expr_, __i );
    }
    inline value_type hessian( int __i, int __j ) const
    {
        return  Value<true/*boost::type_traits::is_fundamental<value_type>::value*/,Expr>::hessian( expr_, __i, __j );
    }
    inline bool deps( int __i ) const
    {
        return expr_.deps( __i ) ;
    }
};

template <class Expr> inline ADExpr< AdFuncExp< ADExpr<Expr> > >
exp ( const ADExpr<Expr>& expr )
{
    typedef AdFuncExp< ADExpr<Expr> > expr_t;
    return ADExpr< expr_t >(  expr_t( expr ) );
}

template <class T, int Nvar, int Order, int Var> inline ADExpr< AdFuncExp< ADType<T, Nvar, Order, Var> > >
exp ( const ADType<T, Nvar, Order, Var>& x )
{
    typedef AdFuncExp< ADType<T, Nvar, Order, Var> > expr_t;
    return ADExpr< expr_t >(  expr_t( x ) );
}
/****************************************************************************
 * Logarithm function
 ****************************************************************************/
/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-02-14

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

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
   \file adfunctions.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-14
 */
template <class Expr> class AdFuncLog
{
public:
    enum { nvar = Expr::nvar };
    typedef typename Expr::value_type value_type;
protected:
    AdFuncLog () {}

    Expr expr_;
public:

    AdFuncLog ( const Expr & expr ) : expr_( expr )
    {
        ;
    }


    template<int isFundamental, typename Expr_>
    struct Value
    {

        typedef typename Expr_::value_type value_type;
        static value_type signum( value_type a )
        {
            if ( abs( a ) < std::numeric_limits<value_type>::epsilon() )
                return value_type( 0 );

            return a/abs( a );
        }
        static value_type value( Expr_ const& expr_ )
        {
            return log( expr_.value() );
        }
        static value_type grad( Expr_ const& expr_, int __i )
        {
            return expr_.grad( __i )/expr_.value();
        }
        static value_type hessian( Expr_ const& expr_, int __i, int __j )
        {
            return expr_.hessian( __i,__j )/expr_.value() - expr_.grad( __i )*expr_.grad( __j )/( pow( expr_.value(),2. ) );
        }
    };
    template<typename Expr_>
    struct Value<true, Expr_>
    {
        typedef typename Expr_::value_type value_type;
        static value_type signum( value_type a )
        {
            if ( std::abs( a ) < std::numeric_limits<value_type>::epsilon() )
                return value_type( 0 );

            return a/std::abs( a );
        }
        static value_type value( Expr_ const& expr_ )
        {
            return std::log( expr_.value() );
        }
        static value_type grad( Expr_ const& expr_, int __i )
        {
            return expr_.grad( __i )/expr_.value();
        }
        static value_type hessian( Expr_ const& expr_, int __i, int __j )
        {
            return expr_.hessian( __i,__j )/expr_.value() - expr_.grad( __i )*expr_.grad( __j )/( std::pow( expr_.value(),2. ) );
        }
    };

    inline value_type value() const
    {
        return Value<true/*boost::type_traits::is_fundamental<value_type>::value*/,Expr>::value( expr_ );
    }
    inline value_type grad( int __i ) const
    {
        return  Value<true/*boost::type_traits::is_fundamental<value_type>::value*/,Expr>::grad( expr_, __i );
    }
    inline value_type hessian( int __i, int __j ) const
    {
        return  Value<true/*boost::type_traits::is_fundamental<value_type>::value*/,Expr>::hessian( expr_, __i, __j );
    }
    inline bool deps( int __i ) const
    {
        return expr_.deps( __i ) ;
    }
};

template <class Expr> inline ADExpr< AdFuncLog< ADExpr<Expr> > >
log ( const ADExpr<Expr>& expr )
{
    typedef AdFuncLog< ADExpr<Expr> > expr_t;
    return ADExpr< expr_t >(  expr_t( expr ) );
}

template <class T, int Nvar, int Order, int Var> inline ADExpr< AdFuncLog< ADType<T, Nvar, Order, Var> > >
log ( const ADType<T, Nvar, Order, Var>& x )
{
    typedef AdFuncLog< ADType<T, Nvar, Order, Var> > expr_t;
    return ADExpr< expr_t >(  expr_t( x ) );
}
/****************************************************************************
 * Log10 function
 ****************************************************************************/
/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-02-14

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

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
   \file adfunctions.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-14
 */
template <class Expr> class AdFuncLog10
{
public:
    enum { nvar = Expr::nvar };
    typedef typename Expr::value_type value_type;
protected:
    AdFuncLog10 () {}

    Expr expr_;
public:

    AdFuncLog10 ( const Expr & expr ) : expr_( expr )
    {
        ;
    }


    template<int isFundamental, typename Expr_>
    struct Value
    {

        typedef typename Expr_::value_type value_type;
        static value_type signum( value_type a )
        {
            if ( abs( a ) < std::numeric_limits<value_type>::epsilon() )
                return value_type( 0 );

            return a/abs( a );
        }
        static value_type value( Expr_ const& expr_ )
        {
            return log10( expr_.value() );
        }
        static value_type grad( Expr_ const& expr_, int __i )
        {
            return expr_.grad( __i )/( log( value_type( 10 ) )*expr_.value() );
        }
        static value_type hessian( Expr_ const& expr_, int __i, int __j )
        {
            return expr_.hessian( __i,__j )/( log( value_type( 10 ) )*expr_.value() ) -  expr_.grad( __i )*expr_.grad( __j )/( log( value_type( 10 ) )*pow( expr_.value(),2. ) );
        }
    };
    template<typename Expr_>
    struct Value<true, Expr_>
    {
        typedef typename Expr_::value_type value_type;
        static value_type signum( value_type a )
        {
            if ( std::abs( a ) < std::numeric_limits<value_type>::epsilon() )
                return value_type( 0 );

            return a/std::abs( a );
        }
        static value_type value( Expr_ const& expr_ )
        {
            return std::log10( expr_.value() );
        }
        static value_type grad( Expr_ const& expr_, int __i )
        {
            return expr_.grad( __i )/( std::log( value_type( 10 ) )*expr_.value() );
        }
        static value_type hessian( Expr_ const& expr_, int __i, int __j )
        {
            return expr_.hessian( __i,__j )/( std::log( value_type( 10 ) )*expr_.value() ) -  expr_.grad( __i )*expr_.grad( __j )/( std::log( value_type( 10 ) )*std::pow( expr_.value(),2. ) );
        }
    };

    inline value_type value() const
    {
        return Value<true/*boost::type_traits::is_fundamental<value_type>::value*/,Expr>::value( expr_ );
    }
    inline value_type grad( int __i ) const
    {
        return  Value<true/*boost::type_traits::is_fundamental<value_type>::value*/,Expr>::grad( expr_, __i );
    }
    inline value_type hessian( int __i, int __j ) const
    {
        return  Value<true/*boost::type_traits::is_fundamental<value_type>::value*/,Expr>::hessian( expr_, __i, __j );
    }
    inline bool deps( int __i ) const
    {
        return expr_.deps( __i ) ;
    }
};

template <class Expr> inline ADExpr< AdFuncLog10< ADExpr<Expr> > >
log10 ( const ADExpr<Expr>& expr )
{
    typedef AdFuncLog10< ADExpr<Expr> > expr_t;
    return ADExpr< expr_t >(  expr_t( expr ) );
}

template <class T, int Nvar, int Order, int Var> inline ADExpr< AdFuncLog10< ADType<T, Nvar, Order, Var> > >
log10 ( const ADType<T, Nvar, Order, Var>& x )
{
    typedef AdFuncLog10< ADType<T, Nvar, Order, Var> > expr_t;
    return ADExpr< expr_t >(  expr_t( x ) );
}
/****************************************************************************
 * Absolute function
 ****************************************************************************/
/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-02-14

  Copyright (C) 2008 Université Joseph Fourier (Grenoble I)

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
   \file adfunctions.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-14
 */
template <class Expr> class AdFuncAbs
{
public:
    enum { nvar = Expr::nvar };
    typedef typename Expr::value_type value_type;
protected:
    AdFuncAbs () {}

    Expr expr_;
public:

    AdFuncAbs ( const Expr & expr ) : expr_( expr )
    {
        ;
    }


    template<int isFundamental, typename Expr_>
    struct Value
    {

        typedef typename Expr_::value_type value_type;
        static value_type signum( value_type a )
        {
            if ( abs( a ) < std::numeric_limits<value_type>::epsilon() )
                return value_type( 0 );

            return a/abs( a );
        }
        static value_type value( Expr_ const& expr_ )
        {
            return abs( expr_.value() );
        }
        static value_type grad( Expr_ const& expr_, int __i )
        {
            return signum( expr_.value() ) * expr_.grad( __i );
        }
        static value_type hessian( Expr_ const& expr_, int __i, int __j )
        {
            return signum( expr_.value() ) * expr_.hessian( __i,__j );
        }
    };
    template<typename Expr_>
    struct Value<true, Expr_>
    {
        typedef typename Expr_::value_type value_type;
        static value_type signum( value_type a )
        {
            if ( std::abs( a ) < std::numeric_limits<value_type>::epsilon() )
                return value_type( 0 );

            return a/std::abs( a );
        }
        static value_type value( Expr_ const& expr_ )
        {
            return std::abs( expr_.value() );
        }
        static value_type grad( Expr_ const& expr_, int __i )
        {
            return signum( expr_.value() ) * expr_.grad( __i );
        }
        static value_type hessian( Expr_ const& expr_, int __i, int __j )
        {
            return signum( expr_.value() ) * expr_.hessian( __i,__j );
        }
    };

    inline value_type value() const
    {
        return Value<true/*boost::type_traits::is_fundamental<value_type>::value*/,Expr>::value( expr_ );
    }
    inline value_type grad( int __i ) const
    {
        return  Value<true/*boost::type_traits::is_fundamental<value_type>::value*/,Expr>::grad( expr_, __i );
    }
    inline value_type hessian( int __i, int __j ) const
    {
        return  Value<true/*boost::type_traits::is_fundamental<value_type>::value*/,Expr>::hessian( expr_, __i, __j );
    }
    inline bool deps( int __i ) const
    {
        return expr_.deps( __i ) ;
    }
};

template <class Expr> inline ADExpr< AdFuncAbs< ADExpr<Expr> > >
abs ( const ADExpr<Expr>& expr )
{
    typedef AdFuncAbs< ADExpr<Expr> > expr_t;
    return ADExpr< expr_t >(  expr_t( expr ) );
}

template <class T, int Nvar, int Order, int Var> inline ADExpr< AdFuncAbs< ADType<T, Nvar, Order, Var> > >
abs ( const ADType<T, Nvar, Order, Var>& x )
{
    typedef AdFuncAbs< ADType<T, Nvar, Order, Var> > expr_t;
    return ADExpr< expr_t >(  expr_t( x ) );
}
}

#endif
