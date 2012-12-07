/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2008-02-14

  Copyright (C) 2008-2011 Université Joseph Fourier (Grenoble I)

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
   \file adoperators.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2008-02-14
 */
#ifndef __ADOperators_H
#define __ADOperators_H 1

#include <feel/feelopt/adexpr.hpp>


namespace Feel
{
//------------------------------- Ad binary operators ------------------------------------------


//------------------------------- Ad addition operators ------------------------------------------
template <class L, class R>
class ADBinaryAdd
{
public:
    enum { nvar = L::nvar };
    typedef typename L::value_type value_type_L;

    typedef typename R::value_type value_type_R;

    typedef typename SNumericalTraits<value_type_L,value_type_R>::promote value_type;

protected:
    ADBinaryAdd() {}

    const L& __left;
    const R& __right;

public:
    ADBinaryAdd( const L& left, const R& rigth ) : __left( left ), __right( rigth )
    {
        ;
    }
    ~ADBinaryAdd()
    {
        ;
    }


    value_type value() const
    {
        return __left.value() + __right.value();
    }
    value_type grad( int __i ) const
    {
        return __left.grad( __i ) + __right.grad( __i );
    }
    value_type hessian( int __i, int __j ) const
    {
        return __left.hessian( __i, __j ) + __right.hessian( __i, __j );
    }


    bool deps( int i ) const
    {
        return __left.deps( i ) || __right.deps( i );
    }
};

#define AD_BIN_ADD_CST(TYPE)                                                          \
template <class L>                                                                     \
class ADBinaryAdd<L, ADCst<TYPE> >                                                   \
{                                                                                      \
public:                                                                                \
   enum { nvar = L::nvar };                                                            \
                                                                                       \
   typedef typename L::value_type value_type_L;                                        \
   typedef TYPE value_type_R;                                                          \
   typedef typename SNumericalTraits<value_type_L,value_type_R>::promote value_type;   \
   typedef ADCst<TYPE> R;                                                             \
                                                                                       \
protected:                                                                             \
   ADBinaryAdd() {}                                                                   \
                                                                                       \
   const L& __left; const  R __right;                                                  \
                                                                                       \
public:                                                                                \
   ADBinaryAdd(const L& left, const R& rigth) : __left(left), __right(rigth) {;}      \
   ~ADBinaryAdd() {;}                                                                 \
                                                                                       \
                                                                                       \
   value_type value() const {return __left.value() + __right.value();}                 \
   value_type grad( int __i ) const {return __left.grad( __i );}                       \
   value_type hessian( int __i, int __j ) const {return __left.hessian(__i, __j);}     \
                                                                                       \
  bool deps( int i ) const { return __left.deps( i ); }                                \
};                                                                                     \
                                                                                       \
                                                                                       \
template <class R>                                                                     \
class ADBinaryAdd< ADCst<TYPE>, R>                                                   \
{                                                                                      \
public:                                                                                \
   enum { nvar = R::nvar };                                                            \
                                                                                       \
   typedef TYPE value_type_L;                                                          \
   typedef typename R::value_type value_type_R;                                        \
   typedef typename SNumericalTraits<value_type_L,value_type_R>::promote value_type;   \
   typedef ADCst<TYPE> L;                                                             \
                                                                                       \
protected:                                                                             \
   ADBinaryAdd() {}                                                                   \
                                                                                       \
   const L __left; const R& __right;                                                   \
                                                                                       \
public:                                                                                \
   ADBinaryAdd(const L& left, const R& rigth) : __left(left), __right(rigth) {;}      \
   ~ADBinaryAdd() {;}                                                                 \
                                                                                       \
                                                                                       \
   value_type value() const {return __left.value() + __right.value();}                 \
   value_type grad( int __i ) const {return __right.grad( __i );}                      \
   value_type hessian( int __i, int __j ) const {return __right.hessian( __i, __j);}   \
                                                                                       \
  bool deps( int i ) const { return __right.deps( i ); }                               \
};

AD_BIN_ADD_CST( int )
AD_BIN_ADD_CST( long int )
AD_BIN_ADD_CST( float )
AD_BIN_ADD_CST( double )
AD_BIN_ADD_CST( long double )

#undef AD_BIN_ADD_CST

//------------------------------- AD substraction operators ------------------------------------------
template <class L, class R>
class ADBinarySubtract
{
public:
    enum { nvar = L::nvar };
    typedef typename L::value_type value_type_L;
    typedef typename R::value_type value_type_R;
    typedef typename SNumericalTraits<value_type_L,value_type_R>::promote value_type;

protected:
    ADBinarySubtract() {}

    const L& __left;
    const R& __right;

public:
    ADBinarySubtract( const L& left, const R& rigth ) : __left( left ), __right( rigth )
    {
        ;
    }
    ~ADBinarySubtract()
    {
        ;
    }


    value_type value() const
    {
        return __left.value() - __right.value();
    }
    value_type grad( int __i ) const
    {
        return __left.grad( __i ) - __right.grad( __i );
    }
    value_type hessian( int __i, int __j ) const
    {
        return __left.hessian( __i, __j ) - __right.hessian( __i, __j );
    }

    bool deps( int i ) const
    {
        return __left.deps( i ) || __right.deps( i );
    }
};

#define AD_BIN_SUB_CST(TYPE)                                                             \
template <class L>                                                                        \
class ADBinarySubtract<L, ADCst<TYPE> >                                                    \
{                                                                                         \
public:                                                                                   \
                                                                                          \
   enum { nvar = L::nvar };                                                               \
                                                                                          \
   typedef typename L::value_type value_type_L;                                           \
   typedef TYPE value_type_R;                                                             \
   typedef typename SNumericalTraits<value_type_L,value_type_R>::promote value_type;      \
   typedef ADCst<TYPE> R;                                                                \
                                                                                          \
                                                                                          \
protected:                                                                                \
   ADBinarySubtract() {}                                                                    \
                                                                                          \
   const L& __left; const R __right;                                                      \
                                                                                          \
public:                                                                                   \
   ADBinarySubtract(const L& left, const R & rigth) : __left(left), __right(rigth) {;}      \
   ~ADBinarySubtract() {;}                                                                  \
                                                                                          \
                                                                                          \
   value_type value() const {return __left.value() - __right.value();}                    \
   value_type grad( int __i) const {return __left.grad(__i);}                             \
   value_type hessian( int __i, int __j ) const {return __left.hessian( __i, __j);}       \
  bool deps( int i ) const { return __left.deps( i ); }                                   \
};                                                                                        \
                                                                                          \
                                                                                          \
template <class R>                                                                        \
class ADBinarySubtract< ADCst<TYPE>, R>                                                    \
{                                                                                         \
public:                                                                                   \
   enum { nvar = R::nvar };                                                               \
                                                                                          \
   typedef TYPE value_type_L;                                                             \
   typedef typename R::value_type value_type_R;                                           \
   typedef typename SNumericalTraits<value_type_L,value_type_R>::promote value_type;      \
   typedef ADCst<TYPE> L;                                                                \
                                                                                          \
protected:                                                                                \
  ADBinarySubtract() {}                                                                     \
                                                                                          \
  const L __left; const R& __right;                                                       \
                                                                                          \
public:                                                                                   \
  ADBinarySubtract(const L& left, const R& rigth) : __left(left), __right(rigth) {;}        \
  ~ADBinarySubtract() {;}                                                                   \
                                                                                          \
                                                                                          \
  value_type value() const {return __left.value() - __right.value();}                     \
  value_type grad(int __i) const {return - __right.grad(__i);}                            \
  value_type hessian( int __i, int __j ) const {return - __right.hessian( __i, __j);}     \
  bool deps( int i ) const { return __right.deps( i ); }                                  \
};

AD_BIN_SUB_CST( int )
AD_BIN_SUB_CST( long int )
AD_BIN_SUB_CST( float )
AD_BIN_SUB_CST( double )
AD_BIN_SUB_CST( long double )

#undef AD_BIN_SUB_CST

//------------------------------- AD multiplication operators ------------------------------------------
template <class L, class R>
class ADBinaryMultiply
{
public:
    enum { nvar = R::nvar };
    typedef typename L::value_type value_type_L;
    typedef typename R::value_type value_type_R;
    typedef typename SNumericalTraits<value_type_L,value_type_R>::promote value_type;

protected:
    ADBinaryMultiply() {}

    const L& __left;
    const R& __right;

public:
    ADBinaryMultiply( const L& left, const R& rigth ) : __left( left ), __right( rigth )
    {
        ;
    }
    ~ADBinaryMultiply()
    {
        ;
    }

    value_type value() const
    {
        return __left.value() * __right.value() ;
    }
    value_type grad( int i ) const
    {
        return  __left.grad( i ) * __right.value() + __right.grad( i ) * __left.value();
    }
    value_type hessian( int __i, int __j ) const
    {
        return ( __left.hessian( __i, __j )*__right.value() + __left.grad( __i ) * __right.grad( __j ) +
                 __left.grad( __j ) * __right.grad( __i ) + __left.value()*__right.hessian( __i, __j ) );
    }

    bool deps( int i ) const
    {
        return __left.deps( i ) || __right.deps( i );
    }
};

#define AD_BIN_MUL_CST(TYPE)                                                                \
template <class L>                                                                           \
class ADBinaryMultiply<L, ADCst<TYPE> >                                                    \
{                                                                                            \
public:                                                                                      \
                                                                                             \
   enum { nvar = L::nvar };                                                                  \
                                                                                             \
   typedef typename L::value_type value_type_L;                                              \
   typedef TYPE value_type_R;                                                                \
   typedef typename SNumericalTraits<value_type_L,value_type_R>::promote value_type;         \
   typedef ADCst<TYPE> R;                                                                   \
                                                                                             \
protected:                                                                                   \
   ADBinaryMultiply() {}                                                                    \
                                                                                             \
   const L& __left; const R __right;                                                         \
                                                                                             \
public:                                                                                      \
   ADBinaryMultiply(const L& left, const R& rigth) : __left(left), __right(rigth) {;}       \
   ~ADBinaryMultiply() {;}                                                                  \
                                                                                             \
   value_type value() const {return __left.value() * __right.value() ;}                      \
   value_type grad(int i) const {return  __left.grad(i) * __right.value();}                  \
   value_type hessian(int i, int j) const {return  __left.hessian(i,j) * __right.value();}   \
   bool deps( int i ) const { return __left.deps( i ); }                                     \
};                                                                                           \
                                                                                             \
template <class R>                                                                           \
class ADBinaryMultiply< ADCst<TYPE>, R>                                                    \
{                                                                                            \
public:                                                                                      \
                                                                                             \
   enum { nvar = R::nvar };                                                                  \
                                                                                             \
   typedef TYPE value_type_L;                                                                \
   typedef typename R::value_type value_type_R;                                              \
   typedef typename SNumericalTraits<value_type_L,value_type_R>::promote value_type;         \
   typedef ADCst<TYPE> L;                                                                   \
                                                                                             \
protected:                                                                                   \
   ADBinaryMultiply() {}                                                                    \
                                                                                             \
   const L __left; const R& __right;                                                         \
                                                                                             \
public:                                                                                      \
   ADBinaryMultiply(const L& left, const R& rigth) : __left(left), __right(rigth) {;}       \
   ~ADBinaryMultiply() {;}                                                                  \
                                                                                             \
   value_type value() const {return __left.value() * __right.value() ;}                      \
   value_type grad(int i) const {return  __right.grad(i) * __left.value();}                  \
   value_type hessian(int i, int j) const {return  __right.hessian(i,j) * __left.value();}   \
   bool deps( int i ) const { return __right.deps( i ); }                                    \
};

AD_BIN_MUL_CST( int )
AD_BIN_MUL_CST( long int )
AD_BIN_MUL_CST( float )
AD_BIN_MUL_CST( double )
AD_BIN_MUL_CST( long double )

#undef AD_BIN_MUL_CST

//------------------------------- AD division operators ------------------------------------------
template <class L, class R> class ADBinaryDivide
{
public:
    enum { nvar = R::nvar };
    typedef typename L::value_type value_type_L;
    typedef typename R::value_type value_type_R;
    typedef typename SNumericalTraits<value_type_L,value_type_R>::promote value_type;

protected:
    ADBinaryDivide() {}

    const L& __left;
    const R& __right;

public:
    ADBinaryDivide( const L& left, const R& rigth ) : __left( left ), __right( rigth )
    {
        ;
    }
    ~ADBinaryDivide()
    {
        ;
    }


    value_type value() const
    {
        return __left.value() / __right.value();
    }
    value_type grad( int i ) const
    {
        return  ( __left.grad( i ) / __right.value() ) - __right.grad( i ) * __left.value() / ( __right.value() * __right.value() ) ;
    }

    value_type hessian( int __i, int __j ) const
    {
        return -( - __left.hessian( __i,__j )*__right.value()*__right.value() +
                  __left.grad( __j )* __right.grad( __i ) * __right.value() +
                  __left.grad( __i )* __right.grad( __j ) * __right.value() -
                  value_type( 2.0 ) * __left.value()* __right.grad( __i ) * __right.grad( __j ) +
                  __left.value() * __right.value() * __right.hessian( __i, __j ) ) / ( __right.value()*__right.value()*__right.value() );
    };

    bool deps( int i ) const
    {
        return __left.deps( i ) || __right.deps( i );
    }
};

#define AD_BIN_DIV_CST(TYPE)                                                                                           \
template <class L>                                                                                                      \
class ADBinaryDivide<L, ADCst<TYPE> >                                                                                 \
{                                                                                                                       \
 public:                                                                                                                \
                                                                                                                        \
   enum { nvar = L::nvar };                                                                                             \
                                                                                                                        \
   typedef typename L::value_type value_type_L;                                                                         \
   typedef TYPE value_type_R;                                                                                           \
   typedef typename SNumericalTraits<value_type_L,value_type_R>::promote value_type;                                    \
   typedef ADCst<TYPE> R;                                                                                              \
                                                                                                                        \
protected:                                                                                                              \
  ADBinaryDivide() {}                                                                                                  \
                                                                                                                        \
  const L& __left; const R __right;                                                                                     \
                                                                                                                        \
 public:                                                                                                                \
  ADBinaryDivide(const L& left, const R& rigth) : __left(left), __right(rigth) {;}                                     \
  ~ADBinaryDivide() {;}                                                                                                \
                                                                                                                        \
                                                                                                                        \
  value_type value() const {return __left.value() / __right.value();}                                                   \
  value_type grad(int i) const {return  __left.grad(i) / __right.value();}                                              \
  value_type hessian(int i, int j) const {return  __left.hessian(i,j) / __right.value();}                               \
  bool deps( int i ) const { return __left.deps( i ); }                                                                 \
};                                                                                                                      \
                                                                                                                        \
template <class R>                                                                                                      \
class ADBinaryDivide< ADCst<TYPE>, R>                                                                                 \
{                                                                                                                       \
public:                                                                                                                 \
                                                                                                                        \
   enum { nvar = R::nvar };                                                                                             \
                                                                                                                        \
   typedef TYPE value_type_L;                                                                                           \
   typedef typename R::value_type value_type_R;                                                                         \
   typedef typename SNumericalTraits<value_type_L,value_type_R>::promote value_type;                                    \
   typedef ADCst<TYPE> L;                                                                                              \
                                                                                                                        \
protected:                                                                                                              \
  ADBinaryDivide() {}                                                                                                  \
                                                                                                                        \
  const L __left; const R& __right;                                                                                     \
                                                                                                                        \
public:                                                                                                                 \
  ADBinaryDivide(L left, const R& right) : __left(left), __right(right) {;}                                            \
  ~ADBinaryDivide() {;}                                                                                                \
                                                                                                                        \
  value_type value() const {return __left.value() / __right.value();}                                                   \
  value_type grad(int i) const {return  (- __right.grad(i) * __left.value() ) / (__right.value() * __right.value()) ;}  \
  value_type hessian(int __i, int __j) const                                                                            \
  { return __left.value() *                                                                                             \
       ( 2 *__right.grad(__i)*__right.grad(__j) - __right.value()*__right.hessian(__i,__j) )/                           \
       std::pow(__right.value(),3.);                                                                                    \
  };                                                                                                                    \
  bool deps( int i ) const { return __right.deps( i ); }                                                                \
};

AD_BIN_DIV_CST( int )
AD_BIN_DIV_CST( long int )
AD_BIN_DIV_CST( float )
AD_BIN_DIV_CST( double )
AD_BIN_DIV_CST( long double )

#undef AD_BIN_DIV_CST


//------------------------------- AD pow function ------------------------------------------
template <class L, class R>
class ADBinaryPow
{
public:
    enum { nvar = R::nvar };
    typedef typename L::value_type value_type_L;
    typedef typename R::value_type value_type_R;

    typedef typename SNumericalTraits<value_type_L,value_type_R>::promote value_type;

protected:
    ADBinaryPow() {}

    const L& __left;
    const R& __right;

public:
    ADBinaryPow( const L& left, const R& rigth ) : __left( left ), __right( rigth )
    {
        ;
    }
    ~ADBinaryPow()
    {
        ;
    }


    value_type value() const
    {
        return std::pow( __left.value(), __right.value() );
    }
    value_type grad( int i ) const
    {
        return  ( __right.grad( i )*std::log( __left.value() )+__right.value()*__left.grad( i )/__left.value() )
                *std::pow( __left.value(), __right.value() );
    }

    bool deps( int i ) const
    {
        return __left.deps( i ) || __right.deps( i );
    }
};


template <class L>
class ADBinaryPow<L, ADCst<typename L::value_type> >
{
public:
    enum { nvar = L::nvar };
    typedef typename L::value_type value_type;
    typedef ADCst<value_type> R;

protected:
    ADBinaryPow() {}

    const L& __left;
    const  R __right;

public:
    ADBinaryPow( const L& left, const R& rigth ) : __left( left ), __right( rigth )
    {
        ;
    }
    ~ADBinaryPow()
    {
        ;
    }


    value_type value() const
    {
        return std::pow( __left.value(),__right.value() ) ;
    }
    value_type grad( int i ) const
    {
        return  ( __right.value()*__left.grad( i )/__left.value() )*std::pow( __left.value(), __right.value() );
    }
    value_type hessian( int i, int j ) const
    {
        //return  std::pow(__left.value(),__right.value())*std::pow(__right.value(),2)*__left.grad(i)*__left.grad(j)/std::pow(__left.value(),2) + std::pow(__left.value(),__right.value())*__right.value()*__left.hessian(i,j)/__left.value() - std::pow(__left.value(),__right.value())*__right.value()*
    }
    bool deps( int i ) const
    {
        return __left.deps( i );
    }
};


template <class R>
class ADBinaryPow< ADCst<typename R::value_type>, R>
{
public:
    enum { nvar = R::nvar };
    typedef typename R::value_type value_type;
    typedef ADCst<value_type> L;

protected:
    ADBinaryPow() {}

    const L __left;
    const R& __right;

public:
    ADBinaryPow( const L& left, const R& rigth ) : __left( left ), __right( rigth )
    {
        ;
    }
    ~ADBinaryPow()
    {
        ;
    }

    const value_type value() const
    {
        return std::pow( __left.value(),__right.value() );
    }
    value_type grad( int i ) const
    {
        return ( __right.grad( i )*std::log( __left.value() ) )*std::pow( __left.value(), __right.value() );
    }

    value_type hessian( int i, int j ) const
    {
        return std::pow( __left.value(),__right.value()-2.0 )*__right.value()*__right.value()*__left.grad( i )*__left.grad( j ) +
               std::pow( __left.value(),__right.value()-1.0 )*__right.value()*__left.hessian( i,j ) -
               std::pow( __left.value(),__right.value()-2.0 )*__right.value()*__left.grad( i )*__left.grad( j );
    }
    bool deps( int i ) const
    {
        return __right.deps( i );
    }
};

template <class L>
class ADBinaryPow< L , ADCst<int> >
{
public:
    enum { nvar = L::nvar };

    typedef typename L::value_type value_type;
    typedef ADCst<int> R;

protected:
    ADBinaryPow() {}

    const L& __left;
    const R __right;

public:
    ADBinaryPow( const L& left, const R& rigth ) : __left( left ), __right( rigth )
    {
        ;
    }
    ~ADBinaryPow()
    {
        ;
    }

    template<int isFundamental, typename L_, typename R_>
    struct Value
    {
        typedef typename L_::value_type value_type;
        static value_type value( L_ const& __left, R_ const& __right )
        {
            return pow( __left.value(),__right.value() );
        }
        static value_type grad( L_ const& __left, R_ const& __right, int __i )
        {
            return value_type( __right.value() )*pow( __left.value(), __right.value()-1 );
        }
        static value_type hessian( L_ const& __left, R_ const& __right, int i, int j )
        {
            return pow( __left.value(),__right.value()-2.0 )*__right.value()*__right.value()*__left.grad( i )*__left.grad( j ) +
                   pow( __left.value(),__right.value()-1.0 )*__right.value()*__left.hessian( i,j ) -
                   pow( __left.value(),__right.value()-2.0 )*__right.value()*__left.grad( i )*__left.grad( j );
        }
    };
    template<typename L_, typename R_>
    struct Value<1, L_, R_>
    {
        typedef typename L_::value_type value_type;
        static value_type value( L_ const& __left, R_ const& __right )
        {
            return std::pow( __left.value(),__right.value() );
        }
        static value_type grad( L_ const& __left, R_ const& __right, int __i )
        {
            return __right.value()*std::pow( __left.value(), __right.value()-1 );
        }
        static value_type hessian( L_ const& __left, R_ const& __right, int i, int j )
        {
            return std::pow( __left.value(),__right.value()-2.0 )*__right.value()*__right.value()*__left.grad( i )*__left.grad( j ) +
                   std::pow( __left.value(),__right.value()-1.0 )*__right.value()*__left.hessian( i,j ) -
                   std::pow( __left.value(),__right.value()-2.0 )*__right.value()*__left.grad( i )*__left.grad( j );
        }
    };
    value_type value() const
    {
        return Value<true/**/,L,R>::value( __left, __right );
    }
    value_type grad( int i ) const
    {
        return Value<true/**/,L,R>::grad( __left, __right, i );
    }
    value_type hessian( int i, int j ) const
    {
        return Value<true/**/,L,R>::hessian( __left, __right, i, j );;
    }
    bool deps( int i ) const
    {
        return __left.deps( i );
    }
};


//------------------------------- AD Greater function ------------------------------------------
#define AD_BIN_REL_CST(Rel, OP,TYPE)                                                                                   \
template <class L>                                                                                                      \
class ADBinary ##Rel <L, ADCst<TYPE> >                                                                                 \
{                                                                                                                       \
 public:                                                                                                                \
                                                                                                                        \
   enum { nvar = L::nvar };                                                                                             \
                                                                                                                        \
   typedef typename L::value_type value_type_L;                                                                         \
   typedef TYPE value_type_R;                                                                                           \
   typedef typename SNumericalTraits<value_type_L,value_type_R>::promote value_type;                                    \
   typedef ADCst<TYPE> R;                                                                                              \
                                                                                                                        \
protected:                                                                                                              \
  ADBinary ##Rel () {}                                                                                                  \
                                                                                                                        \
  const L& __left; const R __right;                                                                                     \
                                                                                                                        \
 public:                                                                                                                \
  ADBinary ##Rel (const L& left, const R& rigth) : __left(left), __right(rigth) {;}                                     \
  ~ADBinary ##Rel () {;}                                                                                                \
                                                                                                                        \
                                                                                                                        \
  value_type value() const { return __left.value() OP __right.value(); }                                                \
  value_type grad(int i) const  { return __left.value() OP __right.value(); }                                           \
  value_type hessian(int i, int j) const { return __left.value() OP __right.value(); }                                   \
  bool deps( int i ) const { return __left.deps( i ); }                                                                 \
};                                                                                                                      \
                                                                                                                        \
template <class R>                                                                                                      \
class ADBinary ##Rel < ADCst<TYPE>, R>                                                                                 \
{                                                                                                                       \
public:                                                                                                                 \
                                                                                                                        \
   enum { nvar = R::nvar };                                                                                             \
                                                                                                                        \
   typedef TYPE value_type_L;                                                                                           \
   typedef typename R::value_type value_type_R;                                                                         \
   typedef typename SNumericalTraits<value_type_L,value_type_R>::promote value_type;                                    \
   typedef ADCst<TYPE> L;                                                                                              \
                                                                                                                        \
protected:                                                                                                              \
  ADBinary ##Rel () {}                                                                                                  \
                                                                                                                        \
  const L __left; const R& __right;                                                                                     \
                                                                                                                        \
public:                                                                                                                 \
  ADBinary ##Rel (L left, const R& right) : __left(left), __right(right) {;}                                            \
  ~ADBinary ##Rel () {;}                                                                                                \
                                                                                                                        \
  value_type value() const { return __left.value() OP __right.value(); }                                                \
  value_type grad(int i) const  { return __left.value() OP __right.value(); }                                           \
  value_type hessian(int i, int j) const { return __left.value() OP __right.value(); }                                  \
  bool deps( int i ) const { return __right.deps( i ); }                                                                \
};

#define AD_BIN_REL( Rel, OP )                                                            \
template <class L, class R>                                                               \
class ADBinary ##Rel                                                                     \
{                                                                                         \
public:                                                                                   \
  enum { nvar = R::nvar };                                                                \
  typedef typename L::value_type value_type_L;                                            \
  typedef typename R::value_type value_type_R;                                            \
                                                                                          \
  typedef typename SNumericalTraits<value_type_L,value_type_R>::promote value_type;       \
                                                                                          \
protected:                                                                                \
  ADBinary ##Rel () {}                                                                   \
                                                                                          \
  const L& __left; const R& __right;                                                      \
                                                                                          \
public:                                                                                   \
  ADBinary ##Rel (const L& left, const R& rigth) : __left(left), __right(rigth) {;}      \
  ~ADBinary ##Rel () {;}                                                                 \
                                                                                          \
                                                                                          \
  value_type value() const {return __left.value() OP __right.value();}                    \
  value_type grad(int i) const  {  return  __left.value() OP __right.value(); }           \
  value_type hessian(int i, int j) const { return __left.value() OP __right.value(); }    \
                                                                                          \
  bool deps( int i ) const { return __left.deps( i ) || __right.deps( i ); }              \
};                                                                                        \
AD_BIN_REL_CST( Rel, OP, int)                                           \
AD_BIN_REL_CST( Rel, OP, long int)                                      \
AD_BIN_REL_CST( Rel, OP, float)                                         \
AD_BIN_REL_CST( Rel, OP, double)                                        \
AD_BIN_REL_CST( Rel, OP, long double)

AD_BIN_REL( Greater, > )
AD_BIN_REL( Less, < )
AD_BIN_REL( GreaterOrEqual, >= )
AD_BIN_REL( LessOrEqual, <= )
AD_BIN_REL( Equal, == )
AD_BIN_REL( NotEqual, != )

#undef AD_BIN_REL
#undef AD_BIN_REL_CST

#include <feel/feelopt/adbinaryoperators.hpp>
}

#endif

