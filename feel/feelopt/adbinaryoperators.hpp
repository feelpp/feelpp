/*
   adbinaryoperators.hpp	ADType expression templates (2 operands)
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

#ifndef AD_BOPS_HPP
#define AD_BOPS_HPP

/****************************************************************************
 * Addition Operators
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
   \file adbinaryoperators.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-14
 */
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4


//------------------------------- AD operators ------------------------------------------

template<class E1, class E2>
inline
ADExpr< ADBinaryAdd< ADExpr<E1>, ADExpr<E2> > >
operator+  ( const ADExpr<E1> &v, const ADExpr<E2> &w )
{
    typedef ADBinaryAdd<ADExpr<E1>, ADExpr<E2> > expr_t;
    return ADExpr<expr_t> ( expr_t ( v , w ) );
}
template<class E, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryAdd<ADExpr<E>,ADType<typename E::value_type, Nvar, order, Var> > >
operator+ ( const ADExpr<E> &e,const ADType<typename E::value_type, Nvar, order, Var>& v )
{
    typedef typename E::value_type A;
    typedef ADBinaryAdd<ADExpr<E>,ADType<typename E::value_type, Nvar, order, Var> > expr_t;
    return ADExpr<expr_t>( expr_t ( e, v ) );
}

template<typename A, int Nvar, int order, int Var1, int Var2>
inline
ADExpr<ADBinaryAdd<ADType<A, Nvar, order, Var1>,ADType<A, Nvar, order, Var2> > >
operator+ ( const ADType<A, Nvar, order, Var1> &e1,const ADType<A, Nvar, order, Var2>& e2 )
{
    typedef ADBinaryAdd<ADType<A, Nvar, order, Var1>,ADType<A, Nvar, order, Var2> > expr_t;
    return ADExpr<expr_t>( expr_t ( e1 , e2 ) );
}

template<class E, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryAdd<ADType<typename E::value_type, Nvar, order, Var>,ADExpr<E> > >
operator+ ( const ADType<typename E::value_type, Nvar, order, Var> &v, const ADExpr<E> &e )
{
    typedef typename E::value_type A;
    typedef ADBinaryAdd<ADType<typename E::value_type, Nvar, order, Var>,ADExpr<E> > expr_t;
    return ADExpr<expr_t> ( expr_t ( v , e ) );
}

template<typename A, typename C, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryAdd<ADCst<C>,ADType<A, Nvar, order, Var> > >
operator+ ( C a, const ADType<A, Nvar, order, Var> &e )
{
    typedef ADBinaryAdd<ADCst<C>,ADType<A, Nvar, order, Var> > expr_t;
    return ADExpr<expr_t> ( expr_t ( ADCst<C>( a ), e  ) );
}

template<typename A, typename C, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryAdd<ADType<A, Nvar, order, Var>,ADCst<C> > >
operator+ ( const ADType<A, Nvar, order, Var> &e, C a )
{
    typedef ADBinaryAdd<ADType<A, Nvar, order, Var>,ADCst<C> > expr_t;
    return ADExpr<expr_t>( expr_t ( e ,ADCst<C>( a ) ) );
}

template<class E, class C>
inline
ADExpr<ADBinaryAdd<ADCst<C>,ADExpr<E> > >
operator+ ( C t, const ADExpr<E> &e )
{
    typedef ADBinaryAdd<ADCst<C>,ADExpr<E> > expr_t;
    return ADExpr<expr_t> ( expr_t ( ADCst<C> ( t ),e ) );
}

template<class E, class C>
inline
ADExpr<ADBinaryAdd<ADExpr<E>,ADCst<C> > >
operator+ ( const ADExpr<E> &e, C t )
{
    typedef ADBinaryAdd<ADExpr<E>,ADCst<C> > expr_t;
    return ADExpr<expr_t>( expr_t ( e, ADCst<C> ( t ) ) );
}

/****************************************************************************
 * Subtraction Operators
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
   \file adbinaryoperators.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-14
 */
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4


//------------------------------- AD operators ------------------------------------------

template<class E1, class E2>
inline
ADExpr< ADBinarySubtract< ADExpr<E1>, ADExpr<E2> > >
operator-  ( const ADExpr<E1> &v, const ADExpr<E2> &w )
{
    typedef ADBinarySubtract<ADExpr<E1>, ADExpr<E2> > expr_t;
    return ADExpr<expr_t> ( expr_t ( v , w ) );
}
template<class E, int Nvar, int order, int Var>
inline
ADExpr<ADBinarySubtract<ADExpr<E>,ADType<typename E::value_type, Nvar, order, Var> > >
operator- ( const ADExpr<E> &e,const ADType<typename E::value_type, Nvar, order, Var>& v )
{
    typedef typename E::value_type A;
    typedef ADBinarySubtract<ADExpr<E>,ADType<typename E::value_type, Nvar, order, Var> > expr_t;
    return ADExpr<expr_t>( expr_t ( e, v ) );
}

template<typename A, int Nvar, int order, int Var1, int Var2>
inline
ADExpr<ADBinarySubtract<ADType<A, Nvar, order, Var1>,ADType<A, Nvar, order, Var2> > >
operator- ( const ADType<A, Nvar, order, Var1> &e1,const ADType<A, Nvar, order, Var2>& e2 )
{
    typedef ADBinarySubtract<ADType<A, Nvar, order, Var1>,ADType<A, Nvar, order, Var2> > expr_t;
    return ADExpr<expr_t>( expr_t ( e1 , e2 ) );
}

template<class E, int Nvar, int order, int Var>
inline
ADExpr<ADBinarySubtract<ADType<typename E::value_type, Nvar, order, Var>,ADExpr<E> > >
operator- ( const ADType<typename E::value_type, Nvar, order, Var> &v, const ADExpr<E> &e )
{
    typedef typename E::value_type A;
    typedef ADBinarySubtract<ADType<typename E::value_type, Nvar, order, Var>,ADExpr<E> > expr_t;
    return ADExpr<expr_t> ( expr_t ( v , e ) );
}

template<typename A, typename C, int Nvar, int order, int Var>
inline
ADExpr<ADBinarySubtract<ADCst<C>,ADType<A, Nvar, order, Var> > >
operator- ( C a, const ADType<A, Nvar, order, Var> &e )
{
    typedef ADBinarySubtract<ADCst<C>,ADType<A, Nvar, order, Var> > expr_t;
    return ADExpr<expr_t> ( expr_t ( ADCst<C>( a ), e  ) );
}

template<typename A, typename C, int Nvar, int order, int Var>
inline
ADExpr<ADBinarySubtract<ADType<A, Nvar, order, Var>,ADCst<C> > >
operator- ( const ADType<A, Nvar, order, Var> &e, C a )
{
    typedef ADBinarySubtract<ADType<A, Nvar, order, Var>,ADCst<C> > expr_t;
    return ADExpr<expr_t>( expr_t ( e ,ADCst<C>( a ) ) );
}

template<class E, class C>
inline
ADExpr<ADBinarySubtract<ADCst<C>,ADExpr<E> > >
operator- ( C t, const ADExpr<E> &e )
{
    typedef ADBinarySubtract<ADCst<C>,ADExpr<E> > expr_t;
    return ADExpr<expr_t> ( expr_t ( ADCst<C> ( t ),e ) );
}

template<class E, class C>
inline
ADExpr<ADBinarySubtract<ADExpr<E>,ADCst<C> > >
operator- ( const ADExpr<E> &e, C t )
{
    typedef ADBinarySubtract<ADExpr<E>,ADCst<C> > expr_t;
    return ADExpr<expr_t>( expr_t ( e, ADCst<C> ( t ) ) );
}

/****************************************************************************
 * Multiplication Operators
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
   \file adbinaryoperators.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-14
 */
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4


//------------------------------- AD operators ------------------------------------------

template<class E1, class E2>
inline
ADExpr< ADBinaryMultiply< ADExpr<E1>, ADExpr<E2> > >
operator*  ( const ADExpr<E1> &v, const ADExpr<E2> &w )
{
    typedef ADBinaryMultiply<ADExpr<E1>, ADExpr<E2> > expr_t;
    return ADExpr<expr_t> ( expr_t ( v , w ) );
}
template<class E, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryMultiply<ADExpr<E>,ADType<typename E::value_type, Nvar, order, Var> > >
operator* ( const ADExpr<E> &e,const ADType<typename E::value_type, Nvar, order, Var>& v )
{
    typedef typename E::value_type A;
    typedef ADBinaryMultiply<ADExpr<E>,ADType<typename E::value_type, Nvar, order, Var> > expr_t;
    return ADExpr<expr_t>( expr_t ( e, v ) );
}

template<typename A, int Nvar, int order, int Var1, int Var2>
inline
ADExpr<ADBinaryMultiply<ADType<A, Nvar, order, Var1>,ADType<A, Nvar, order, Var2> > >
operator* ( const ADType<A, Nvar, order, Var1> &e1,const ADType<A, Nvar, order, Var2>& e2 )
{
    typedef ADBinaryMultiply<ADType<A, Nvar, order, Var1>,ADType<A, Nvar, order, Var2> > expr_t;
    return ADExpr<expr_t>( expr_t ( e1 , e2 ) );
}

template<class E, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryMultiply<ADType<typename E::value_type, Nvar, order, Var>,ADExpr<E> > >
operator* ( const ADType<typename E::value_type, Nvar, order, Var> &v, const ADExpr<E> &e )
{
    typedef typename E::value_type A;
    typedef ADBinaryMultiply<ADType<typename E::value_type, Nvar, order, Var>,ADExpr<E> > expr_t;
    return ADExpr<expr_t> ( expr_t ( v , e ) );
}

template<typename A, typename C, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryMultiply<ADCst<C>,ADType<A, Nvar, order, Var> > >
operator* ( C a, const ADType<A, Nvar, order, Var> &e )
{
    typedef ADBinaryMultiply<ADCst<C>,ADType<A, Nvar, order, Var> > expr_t;
    return ADExpr<expr_t> ( expr_t ( ADCst<C>( a ), e  ) );
}

template<typename A, typename C, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryMultiply<ADType<A, Nvar, order, Var>,ADCst<C> > >
operator* ( const ADType<A, Nvar, order, Var> &e, C a )
{
    typedef ADBinaryMultiply<ADType<A, Nvar, order, Var>,ADCst<C> > expr_t;
    return ADExpr<expr_t>( expr_t ( e ,ADCst<C>( a ) ) );
}

template<class E, class C>
inline
ADExpr<ADBinaryMultiply<ADCst<C>,ADExpr<E> > >
operator* ( C t, const ADExpr<E> &e )
{
    typedef ADBinaryMultiply<ADCst<C>,ADExpr<E> > expr_t;
    return ADExpr<expr_t> ( expr_t ( ADCst<C> ( t ),e ) );
}

template<class E, class C>
inline
ADExpr<ADBinaryMultiply<ADExpr<E>,ADCst<C> > >
operator* ( const ADExpr<E> &e, C t )
{
    typedef ADBinaryMultiply<ADExpr<E>,ADCst<C> > expr_t;
    return ADExpr<expr_t>( expr_t ( e, ADCst<C> ( t ) ) );
}

/****************************************************************************
 * Division Operators
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
   \file adbinaryoperators.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-14
 */
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4


//------------------------------- AD operators ------------------------------------------

template<class E1, class E2>
inline
ADExpr< ADBinaryDivide< ADExpr<E1>, ADExpr<E2> > >
operator/  ( const ADExpr<E1> &v, const ADExpr<E2> &w )
{
    typedef ADBinaryDivide<ADExpr<E1>, ADExpr<E2> > expr_t;
    return ADExpr<expr_t> ( expr_t ( v , w ) );
}
template<class E, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryDivide<ADExpr<E>,ADType<typename E::value_type, Nvar, order, Var> > >
operator/ ( const ADExpr<E> &e,const ADType<typename E::value_type, Nvar, order, Var>& v )
{
    typedef typename E::value_type A;
    typedef ADBinaryDivide<ADExpr<E>,ADType<typename E::value_type, Nvar, order, Var> > expr_t;
    return ADExpr<expr_t>( expr_t ( e, v ) );
}

template<typename A, int Nvar, int order, int Var1, int Var2>
inline
ADExpr<ADBinaryDivide<ADType<A, Nvar, order, Var1>,ADType<A, Nvar, order, Var2> > >
operator/ ( const ADType<A, Nvar, order, Var1> &e1,const ADType<A, Nvar, order, Var2>& e2 )
{
    typedef ADBinaryDivide<ADType<A, Nvar, order, Var1>,ADType<A, Nvar, order, Var2> > expr_t;
    return ADExpr<expr_t>( expr_t ( e1 , e2 ) );
}

template<class E, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryDivide<ADType<typename E::value_type, Nvar, order, Var>,ADExpr<E> > >
operator/ ( const ADType<typename E::value_type, Nvar, order, Var> &v, const ADExpr<E> &e )
{
    typedef typename E::value_type A;
    typedef ADBinaryDivide<ADType<typename E::value_type, Nvar, order, Var>,ADExpr<E> > expr_t;
    return ADExpr<expr_t> ( expr_t ( v , e ) );
}

template<typename A, typename C, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryDivide<ADCst<C>,ADType<A, Nvar, order, Var> > >
operator/ ( C a, const ADType<A, Nvar, order, Var> &e )
{
    typedef ADBinaryDivide<ADCst<C>,ADType<A, Nvar, order, Var> > expr_t;
    return ADExpr<expr_t> ( expr_t ( ADCst<C>( a ), e  ) );
}

template<typename A, typename C, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryDivide<ADType<A, Nvar, order, Var>,ADCst<C> > >
operator/ ( const ADType<A, Nvar, order, Var> &e, C a )
{
    typedef ADBinaryDivide<ADType<A, Nvar, order, Var>,ADCst<C> > expr_t;
    return ADExpr<expr_t>( expr_t ( e ,ADCst<C>( a ) ) );
}

template<class E, class C>
inline
ADExpr<ADBinaryDivide<ADCst<C>,ADExpr<E> > >
operator/ ( C t, const ADExpr<E> &e )
{
    typedef ADBinaryDivide<ADCst<C>,ADExpr<E> > expr_t;
    return ADExpr<expr_t> ( expr_t ( ADCst<C> ( t ),e ) );
}

template<class E, class C>
inline
ADExpr<ADBinaryDivide<ADExpr<E>,ADCst<C> > >
operator/ ( const ADExpr<E> &e, C t )
{
    typedef ADBinaryDivide<ADExpr<E>,ADCst<C> > expr_t;
    return ADExpr<expr_t>( expr_t ( e, ADCst<C> ( t ) ) );
}

/****************************************************************************
 * Pow Operators
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
   \file adbinaryoperators.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-14
 */
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4


//------------------------------- AD operators ------------------------------------------

template<class E1, class E2>
inline
ADExpr< ADBinaryPow< ADExpr<E1>, ADExpr<E2> > >
operator^  ( const ADExpr<E1> &v, const ADExpr<E2> &w )
{
    typedef ADBinaryPow<ADExpr<E1>, ADExpr<E2> > expr_t;
    return ADExpr<expr_t> ( expr_t ( v , w ) );
}
template<class E, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryPow<ADExpr<E>,ADType<typename E::value_type, Nvar, order, Var> > >
operator^ ( const ADExpr<E> &e,const ADType<typename E::value_type, Nvar, order, Var>& v )
{
    typedef typename E::value_type A;
    typedef ADBinaryPow<ADExpr<E>,ADType<typename E::value_type, Nvar, order, Var> > expr_t;
    return ADExpr<expr_t>( expr_t ( e, v ) );
}

template<typename A, int Nvar, int order, int Var1, int Var2>
inline
ADExpr<ADBinaryPow<ADType<A, Nvar, order, Var1>,ADType<A, Nvar, order, Var2> > >
operator^ ( const ADType<A, Nvar, order, Var1> &e1,const ADType<A, Nvar, order, Var2>& e2 )
{
    typedef ADBinaryPow<ADType<A, Nvar, order, Var1>,ADType<A, Nvar, order, Var2> > expr_t;
    return ADExpr<expr_t>( expr_t ( e1 , e2 ) );
}

template<class E, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryPow<ADType<typename E::value_type, Nvar, order, Var>,ADExpr<E> > >
operator^ ( const ADType<typename E::value_type, Nvar, order, Var> &v, const ADExpr<E> &e )
{
    typedef typename E::value_type A;
    typedef ADBinaryPow<ADType<typename E::value_type, Nvar, order, Var>,ADExpr<E> > expr_t;
    return ADExpr<expr_t> ( expr_t ( v , e ) );
}

template<typename A, typename C, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryPow<ADCst<C>,ADType<A, Nvar, order, Var> > >
operator^ ( C a, const ADType<A, Nvar, order, Var> &e )
{
    typedef ADBinaryPow<ADCst<C>,ADType<A, Nvar, order, Var> > expr_t;
    return ADExpr<expr_t> ( expr_t ( ADCst<C>( a ), e  ) );
}

template<typename A, typename C, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryPow<ADType<A, Nvar, order, Var>,ADCst<C> > >
operator^ ( const ADType<A, Nvar, order, Var> &e, C a )
{
    typedef ADBinaryPow<ADType<A, Nvar, order, Var>,ADCst<C> > expr_t;
    return ADExpr<expr_t>( expr_t ( e ,ADCst<C>( a ) ) );
}

template<class E, class C>
inline
ADExpr<ADBinaryPow<ADCst<C>,ADExpr<E> > >
operator^ ( C t, const ADExpr<E> &e )
{
    typedef ADBinaryPow<ADCst<C>,ADExpr<E> > expr_t;
    return ADExpr<expr_t> ( expr_t ( ADCst<C> ( t ),e ) );
}

template<class E, class C>
inline
ADExpr<ADBinaryPow<ADExpr<E>,ADCst<C> > >
operator^ ( const ADExpr<E> &e, C t )
{
    typedef ADBinaryPow<ADExpr<E>,ADCst<C> > expr_t;
    return ADExpr<expr_t>( expr_t ( e, ADCst<C> ( t ) ) );
}

/****************************************************************************
 * Greater-than Operators
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
   \file adbinaryoperators.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-14
 */
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4


//------------------------------- AD operators ------------------------------------------

template<class E1, class E2>
inline
ADExpr< ADBinaryGreater< ADExpr<E1>, ADExpr<E2> > >
operator>  ( const ADExpr<E1> &v, const ADExpr<E2> &w )
{
    typedef ADBinaryGreater<ADExpr<E1>, ADExpr<E2> > expr_t;
    return ADExpr<expr_t> ( expr_t ( v , w ) );
}
template<class E, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryGreater<ADExpr<E>,ADType<typename E::value_type, Nvar, order, Var> > >
operator> ( const ADExpr<E> &e,const ADType<typename E::value_type, Nvar, order, Var>& v )
{
    typedef typename E::value_type A;
    typedef ADBinaryGreater<ADExpr<E>,ADType<typename E::value_type, Nvar, order, Var> > expr_t;
    return ADExpr<expr_t>( expr_t ( e, v ) );
}

template<typename A, int Nvar, int order, int Var1, int Var2>
inline
ADExpr<ADBinaryGreater<ADType<A, Nvar, order, Var1>,ADType<A, Nvar, order, Var2> > >
operator> ( const ADType<A, Nvar, order, Var1> &e1,const ADType<A, Nvar, order, Var2>& e2 )
{
    typedef ADBinaryGreater<ADType<A, Nvar, order, Var1>,ADType<A, Nvar, order, Var2> > expr_t;
    return ADExpr<expr_t>( expr_t ( e1 , e2 ) );
}

template<class E, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryGreater<ADType<typename E::value_type, Nvar, order, Var>,ADExpr<E> > >
operator> ( const ADType<typename E::value_type, Nvar, order, Var> &v, const ADExpr<E> &e )
{
    typedef typename E::value_type A;
    typedef ADBinaryGreater<ADType<typename E::value_type, Nvar, order, Var>,ADExpr<E> > expr_t;
    return ADExpr<expr_t> ( expr_t ( v , e ) );
}

template<typename A, typename C, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryGreater<ADCst<C>,ADType<A, Nvar, order, Var> > >
operator> ( C a, const ADType<A, Nvar, order, Var> &e )
{
    typedef ADBinaryGreater<ADCst<C>,ADType<A, Nvar, order, Var> > expr_t;
    return ADExpr<expr_t> ( expr_t ( ADCst<C>( a ), e  ) );
}

template<typename A, typename C, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryGreater<ADType<A, Nvar, order, Var>,ADCst<C> > >
operator> ( const ADType<A, Nvar, order, Var> &e, C a )
{
    typedef ADBinaryGreater<ADType<A, Nvar, order, Var>,ADCst<C> > expr_t;
    return ADExpr<expr_t>( expr_t ( e ,ADCst<C>( a ) ) );
}

template<class E, class C>
inline
ADExpr<ADBinaryGreater<ADCst<C>,ADExpr<E> > >
operator> ( C t, const ADExpr<E> &e )
{
    typedef ADBinaryGreater<ADCst<C>,ADExpr<E> > expr_t;
    return ADExpr<expr_t> ( expr_t ( ADCst<C> ( t ),e ) );
}

template<class E, class C>
inline
ADExpr<ADBinaryGreater<ADExpr<E>,ADCst<C> > >
operator> ( const ADExpr<E> &e, C t )
{
    typedef ADBinaryGreater<ADExpr<E>,ADCst<C> > expr_t;
    return ADExpr<expr_t>( expr_t ( e, ADCst<C> ( t ) ) );
}

/****************************************************************************
 * Less-than Operators
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
   \file adbinaryoperators.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-14
 */
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4


//------------------------------- AD operators ------------------------------------------

template<class E1, class E2>
inline
ADExpr< ADBinaryLess< ADExpr<E1>, ADExpr<E2> > >
operator<  ( const ADExpr<E1> &v, const ADExpr<E2> &w )
{
    typedef ADBinaryLess<ADExpr<E1>, ADExpr<E2> > expr_t;
    return ADExpr<expr_t> ( expr_t ( v , w ) );
}
template<class E, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryLess<ADExpr<E>,ADType<typename E::value_type, Nvar, order, Var> > >
operator< ( const ADExpr<E> &e,const ADType<typename E::value_type, Nvar, order, Var>& v )
{
    typedef typename E::value_type A;
    typedef ADBinaryLess<ADExpr<E>,ADType<typename E::value_type, Nvar, order, Var> > expr_t;
    return ADExpr<expr_t>( expr_t ( e, v ) );
}

template<typename A, int Nvar, int order, int Var1, int Var2>
inline
ADExpr<ADBinaryLess<ADType<A, Nvar, order, Var1>,ADType<A, Nvar, order, Var2> > >
operator< ( const ADType<A, Nvar, order, Var1> &e1,const ADType<A, Nvar, order, Var2>& e2 )
{
    typedef ADBinaryLess<ADType<A, Nvar, order, Var1>,ADType<A, Nvar, order, Var2> > expr_t;
    return ADExpr<expr_t>( expr_t ( e1 , e2 ) );
}

template<class E, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryLess<ADType<typename E::value_type, Nvar, order, Var>,ADExpr<E> > >
operator< ( const ADType<typename E::value_type, Nvar, order, Var> &v, const ADExpr<E> &e )
{
    typedef typename E::value_type A;
    typedef ADBinaryLess<ADType<typename E::value_type, Nvar, order, Var>,ADExpr<E> > expr_t;
    return ADExpr<expr_t> ( expr_t ( v , e ) );
}

template<typename A, typename C, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryLess<ADCst<C>,ADType<A, Nvar, order, Var> > >
operator< ( C a, const ADType<A, Nvar, order, Var> &e )
{
    typedef ADBinaryLess<ADCst<C>,ADType<A, Nvar, order, Var> > expr_t;
    return ADExpr<expr_t> ( expr_t ( ADCst<C>( a ), e  ) );
}

template<typename A, typename C, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryLess<ADType<A, Nvar, order, Var>,ADCst<C> > >
operator< ( const ADType<A, Nvar, order, Var> &e, C a )
{
    typedef ADBinaryLess<ADType<A, Nvar, order, Var>,ADCst<C> > expr_t;
    return ADExpr<expr_t>( expr_t ( e ,ADCst<C>( a ) ) );
}

template<class E, class C>
inline
ADExpr<ADBinaryLess<ADCst<C>,ADExpr<E> > >
operator< ( C t, const ADExpr<E> &e )
{
    typedef ADBinaryLess<ADCst<C>,ADExpr<E> > expr_t;
    return ADExpr<expr_t> ( expr_t ( ADCst<C> ( t ),e ) );
}

template<class E, class C>
inline
ADExpr<ADBinaryLess<ADExpr<E>,ADCst<C> > >
operator< ( const ADExpr<E> &e, C t )
{
    typedef ADBinaryLess<ADExpr<E>,ADCst<C> > expr_t;
    return ADExpr<expr_t>( expr_t ( e, ADCst<C> ( t ) ) );
}

/****************************************************************************
 * Greater or equal (>=) operators
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
   \file adbinaryoperators.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-14
 */
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4


//------------------------------- AD operators ------------------------------------------

template<class E1, class E2>
inline
ADExpr< ADBinaryGreaterOrEqual< ADExpr<E1>, ADExpr<E2> > >
operator>=  ( const ADExpr<E1> &v, const ADExpr<E2> &w )
{
    typedef ADBinaryGreaterOrEqual<ADExpr<E1>, ADExpr<E2> > expr_t;
    return ADExpr<expr_t> ( expr_t ( v , w ) );
}
template<class E, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryGreaterOrEqual<ADExpr<E>,ADType<typename E::value_type, Nvar, order, Var> > >
operator>= ( const ADExpr<E> &e,const ADType<typename E::value_type, Nvar, order, Var>& v )
{
    typedef typename E::value_type A;
    typedef ADBinaryGreaterOrEqual<ADExpr<E>,ADType<typename E::value_type, Nvar, order, Var> > expr_t;
    return ADExpr<expr_t>( expr_t ( e, v ) );
}

template<typename A, int Nvar, int order, int Var1, int Var2>
inline
ADExpr<ADBinaryGreaterOrEqual<ADType<A, Nvar, order, Var1>,ADType<A, Nvar, order, Var2> > >
operator>= ( const ADType<A, Nvar, order, Var1> &e1,const ADType<A, Nvar, order, Var2>& e2 )
{
    typedef ADBinaryGreaterOrEqual<ADType<A, Nvar, order, Var1>,ADType<A, Nvar, order, Var2> > expr_t;
    return ADExpr<expr_t>( expr_t ( e1 , e2 ) );
}

template<class E, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryGreaterOrEqual<ADType<typename E::value_type, Nvar, order, Var>,ADExpr<E> > >
operator>= ( const ADType<typename E::value_type, Nvar, order, Var> &v, const ADExpr<E> &e )
{
    typedef typename E::value_type A;
    typedef ADBinaryGreaterOrEqual<ADType<typename E::value_type, Nvar, order, Var>,ADExpr<E> > expr_t;
    return ADExpr<expr_t> ( expr_t ( v , e ) );
}

template<typename A, typename C, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryGreaterOrEqual<ADCst<C>,ADType<A, Nvar, order, Var> > >
operator>= ( C a, const ADType<A, Nvar, order, Var> &e )
{
    typedef ADBinaryGreaterOrEqual<ADCst<C>,ADType<A, Nvar, order, Var> > expr_t;
    return ADExpr<expr_t> ( expr_t ( ADCst<C>( a ), e  ) );
}

template<typename A, typename C, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryGreaterOrEqual<ADType<A, Nvar, order, Var>,ADCst<C> > >
operator>= ( const ADType<A, Nvar, order, Var> &e, C a )
{
    typedef ADBinaryGreaterOrEqual<ADType<A, Nvar, order, Var>,ADCst<C> > expr_t;
    return ADExpr<expr_t>( expr_t ( e ,ADCst<C>( a ) ) );
}

template<class E, class C>
inline
ADExpr<ADBinaryGreaterOrEqual<ADCst<C>,ADExpr<E> > >
operator>= ( C t, const ADExpr<E> &e )
{
    typedef ADBinaryGreaterOrEqual<ADCst<C>,ADExpr<E> > expr_t;
    return ADExpr<expr_t> ( expr_t ( ADCst<C> ( t ),e ) );
}

template<class E, class C>
inline
ADExpr<ADBinaryGreaterOrEqual<ADExpr<E>,ADCst<C> > >
operator>= ( const ADExpr<E> &e, C t )
{
    typedef ADBinaryGreaterOrEqual<ADExpr<E>,ADCst<C> > expr_t;
    return ADExpr<expr_t>( expr_t ( e, ADCst<C> ( t ) ) );
}

/****************************************************************************
 * Less or equal (<=) operators
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
   \file adbinaryoperators.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-14
 */
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4


//------------------------------- AD operators ------------------------------------------

template<class E1, class E2>
inline
ADExpr< ADBinaryLessOrEqual< ADExpr<E1>, ADExpr<E2> > >
operator<=  ( const ADExpr<E1> &v, const ADExpr<E2> &w )
{
    typedef ADBinaryLessOrEqual<ADExpr<E1>, ADExpr<E2> > expr_t;
    return ADExpr<expr_t> ( expr_t ( v , w ) );
}
template<class E, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryLessOrEqual<ADExpr<E>,ADType<typename E::value_type, Nvar, order, Var> > >
operator<= ( const ADExpr<E> &e,const ADType<typename E::value_type, Nvar, order, Var>& v )
{
    typedef typename E::value_type A;
    typedef ADBinaryLessOrEqual<ADExpr<E>,ADType<typename E::value_type, Nvar, order, Var> > expr_t;
    return ADExpr<expr_t>( expr_t ( e, v ) );
}

template<typename A, int Nvar, int order, int Var1, int Var2>
inline
ADExpr<ADBinaryLessOrEqual<ADType<A, Nvar, order, Var1>,ADType<A, Nvar, order, Var2> > >
operator<= ( const ADType<A, Nvar, order, Var1> &e1,const ADType<A, Nvar, order, Var2>& e2 )
{
    typedef ADBinaryLessOrEqual<ADType<A, Nvar, order, Var1>,ADType<A, Nvar, order, Var2> > expr_t;
    return ADExpr<expr_t>( expr_t ( e1 , e2 ) );
}

template<class E, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryLessOrEqual<ADType<typename E::value_type, Nvar, order, Var>,ADExpr<E> > >
operator<= ( const ADType<typename E::value_type, Nvar, order, Var> &v, const ADExpr<E> &e )
{
    typedef typename E::value_type A;
    typedef ADBinaryLessOrEqual<ADType<typename E::value_type, Nvar, order, Var>,ADExpr<E> > expr_t;
    return ADExpr<expr_t> ( expr_t ( v , e ) );
}

template<typename A, typename C, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryLessOrEqual<ADCst<C>,ADType<A, Nvar, order, Var> > >
operator<= ( C a, const ADType<A, Nvar, order, Var> &e )
{
    typedef ADBinaryLessOrEqual<ADCst<C>,ADType<A, Nvar, order, Var> > expr_t;
    return ADExpr<expr_t> ( expr_t ( ADCst<C>( a ), e  ) );
}

template<typename A, typename C, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryLessOrEqual<ADType<A, Nvar, order, Var>,ADCst<C> > >
operator<= ( const ADType<A, Nvar, order, Var> &e, C a )
{
    typedef ADBinaryLessOrEqual<ADType<A, Nvar, order, Var>,ADCst<C> > expr_t;
    return ADExpr<expr_t>( expr_t ( e ,ADCst<C>( a ) ) );
}

template<class E, class C>
inline
ADExpr<ADBinaryLessOrEqual<ADCst<C>,ADExpr<E> > >
operator<= ( C t, const ADExpr<E> &e )
{
    typedef ADBinaryLessOrEqual<ADCst<C>,ADExpr<E> > expr_t;
    return ADExpr<expr_t> ( expr_t ( ADCst<C> ( t ),e ) );
}

template<class E, class C>
inline
ADExpr<ADBinaryLessOrEqual<ADExpr<E>,ADCst<C> > >
operator<= ( const ADExpr<E> &e, C t )
{
    typedef ADBinaryLessOrEqual<ADExpr<E>,ADCst<C> > expr_t;
    return ADExpr<expr_t>( expr_t ( e, ADCst<C> ( t ) ) );
}

/****************************************************************************
 * Equality operators
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
   \file adbinaryoperators.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-14
 */
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4


//------------------------------- AD operators ------------------------------------------

template<class E1, class E2>
inline
ADExpr< ADBinaryEqual< ADExpr<E1>, ADExpr<E2> > >
operator==  ( const ADExpr<E1> &v, const ADExpr<E2> &w )
{
    typedef ADBinaryEqual<ADExpr<E1>, ADExpr<E2> > expr_t;
    return ADExpr<expr_t> ( expr_t ( v , w ) );
}
template<class E, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryEqual<ADExpr<E>,ADType<typename E::value_type, Nvar, order, Var> > >
operator== ( const ADExpr<E> &e,const ADType<typename E::value_type, Nvar, order, Var>& v )
{
    typedef typename E::value_type A;
    typedef ADBinaryEqual<ADExpr<E>,ADType<typename E::value_type, Nvar, order, Var> > expr_t;
    return ADExpr<expr_t>( expr_t ( e, v ) );
}

template<typename A, int Nvar, int order, int Var1, int Var2>
inline
ADExpr<ADBinaryEqual<ADType<A, Nvar, order, Var1>,ADType<A, Nvar, order, Var2> > >
operator== ( const ADType<A, Nvar, order, Var1> &e1,const ADType<A, Nvar, order, Var2>& e2 )
{
    typedef ADBinaryEqual<ADType<A, Nvar, order, Var1>,ADType<A, Nvar, order, Var2> > expr_t;
    return ADExpr<expr_t>( expr_t ( e1 , e2 ) );
}

template<class E, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryEqual<ADType<typename E::value_type, Nvar, order, Var>,ADExpr<E> > >
operator== ( const ADType<typename E::value_type, Nvar, order, Var> &v, const ADExpr<E> &e )
{
    typedef typename E::value_type A;
    typedef ADBinaryEqual<ADType<typename E::value_type, Nvar, order, Var>,ADExpr<E> > expr_t;
    return ADExpr<expr_t> ( expr_t ( v , e ) );
}

template<typename A, typename C, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryEqual<ADCst<C>,ADType<A, Nvar, order, Var> > >
operator== ( C a, const ADType<A, Nvar, order, Var> &e )
{
    typedef ADBinaryEqual<ADCst<C>,ADType<A, Nvar, order, Var> > expr_t;
    return ADExpr<expr_t> ( expr_t ( ADCst<C>( a ), e  ) );
}

template<typename A, typename C, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryEqual<ADType<A, Nvar, order, Var>,ADCst<C> > >
operator== ( const ADType<A, Nvar, order, Var> &e, C a )
{
    typedef ADBinaryEqual<ADType<A, Nvar, order, Var>,ADCst<C> > expr_t;
    return ADExpr<expr_t>( expr_t ( e ,ADCst<C>( a ) ) );
}

template<class E, class C>
inline
ADExpr<ADBinaryEqual<ADCst<C>,ADExpr<E> > >
operator== ( C t, const ADExpr<E> &e )
{
    typedef ADBinaryEqual<ADCst<C>,ADExpr<E> > expr_t;
    return ADExpr<expr_t> ( expr_t ( ADCst<C> ( t ),e ) );
}

template<class E, class C>
inline
ADExpr<ADBinaryEqual<ADExpr<E>,ADCst<C> > >
operator== ( const ADExpr<E> &e, C t )
{
    typedef ADBinaryEqual<ADExpr<E>,ADCst<C> > expr_t;
    return ADExpr<expr_t>( expr_t ( e, ADCst<C> ( t ) ) );
}

/****************************************************************************
 * Not-equal operators
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
   \file adbinaryoperators.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-14
 */
// -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4


//------------------------------- AD operators ------------------------------------------

template<class E1, class E2>
inline
ADExpr< ADBinaryNotEqual< ADExpr<E1>, ADExpr<E2> > >
operator!=  ( const ADExpr<E1> &v, const ADExpr<E2> &w )
{
    typedef ADBinaryNotEqual<ADExpr<E1>, ADExpr<E2> > expr_t;
    return ADExpr<expr_t> ( expr_t ( v , w ) );
}
template<class E, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryNotEqual<ADExpr<E>,ADType<typename E::value_type, Nvar, order, Var> > >
operator!= ( const ADExpr<E> &e,const ADType<typename E::value_type, Nvar, order, Var>& v )
{
    typedef typename E::value_type A;
    typedef ADBinaryNotEqual<ADExpr<E>,ADType<typename E::value_type, Nvar, order, Var> > expr_t;
    return ADExpr<expr_t>( expr_t ( e, v ) );
}

template<typename A, int Nvar, int order, int Var1, int Var2>
inline
ADExpr<ADBinaryNotEqual<ADType<A, Nvar, order, Var1>,ADType<A, Nvar, order, Var2> > >
operator!= ( const ADType<A, Nvar, order, Var1> &e1,const ADType<A, Nvar, order, Var2>& e2 )
{
    typedef ADBinaryNotEqual<ADType<A, Nvar, order, Var1>,ADType<A, Nvar, order, Var2> > expr_t;
    return ADExpr<expr_t>( expr_t ( e1 , e2 ) );
}

template<class E, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryNotEqual<ADType<typename E::value_type, Nvar, order, Var>,ADExpr<E> > >
operator!= ( const ADType<typename E::value_type, Nvar, order, Var> &v, const ADExpr<E> &e )
{
    typedef typename E::value_type A;
    typedef ADBinaryNotEqual<ADType<typename E::value_type, Nvar, order, Var>,ADExpr<E> > expr_t;
    return ADExpr<expr_t> ( expr_t ( v , e ) );
}

template<typename A, typename C, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryNotEqual<ADCst<C>,ADType<A, Nvar, order, Var> > >
operator!= ( C a, const ADType<A, Nvar, order, Var> &e )
{
    typedef ADBinaryNotEqual<ADCst<C>,ADType<A, Nvar, order, Var> > expr_t;
    return ADExpr<expr_t> ( expr_t ( ADCst<C>( a ), e  ) );
}

template<typename A, typename C, int Nvar, int order, int Var>
inline
ADExpr<ADBinaryNotEqual<ADType<A, Nvar, order, Var>,ADCst<C> > >
operator!= ( const ADType<A, Nvar, order, Var> &e, C a )
{
    typedef ADBinaryNotEqual<ADType<A, Nvar, order, Var>,ADCst<C> > expr_t;
    return ADExpr<expr_t>( expr_t ( e ,ADCst<C>( a ) ) );
}

template<class E, class C>
inline
ADExpr<ADBinaryNotEqual<ADCst<C>,ADExpr<E> > >
operator!= ( C t, const ADExpr<E> &e )
{
    typedef ADBinaryNotEqual<ADCst<C>,ADExpr<E> > expr_t;
    return ADExpr<expr_t> ( expr_t ( ADCst<C> ( t ),e ) );
}

template<class E, class C>
inline
ADExpr<ADBinaryNotEqual<ADExpr<E>,ADCst<C> > >
operator!= ( const ADExpr<E> &e, C t )
{
    typedef ADBinaryNotEqual<ADExpr<E>,ADCst<C> > expr_t;
    return ADExpr<expr_t>( expr_t ( e, ADCst<C> ( t ) ) );
}


#endif
