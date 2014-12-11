/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2007-05-30

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
   \file form.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2007-05-30
 */
#ifndef __Form_H
#define __Form_H 1

#include <feel/feelcore/parameter.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feelalg/vector.hpp>
#include <feel/feelalg/matrixsparse.hpp>
#include <feel/feelalg/backend.hpp>
#include <feel/feeldiscr/functionspace.hpp>

#include <feel/feelvf/bilinearform.hpp>
#include <feel/feelvf/linearform.hpp>
#include <feel/feelvf/integrate.hpp>

namespace Feel
{
template<typename T> class Vector;

//
// free form functions
//
template<typename X1, typename X2>
inline
vf::detail::BilinearForm<X1, X2>
form( boost::shared_ptr<X1> const& __X1,
      boost::shared_ptr<X2> const& __X2,
      boost::shared_ptr<MatrixSparse<double> > __M,
      size_type rowstart = 0,
      size_type colstart = 0,
      bool init = false,
      bool do_threshold = false,
      typename X1::value_type threshold = type_traits<double>::epsilon(),
      size_type pattern  = Pattern::COUPLED )
{
    return vf::detail::BilinearForm<X1, X2>( __X1, __X2, __M, rowstart, colstart, init, do_threshold, threshold, pattern );
}

template<typename X1, typename RepType>
inline
vf::detail::LinearForm<X1, RepType, RepType>
form( boost::shared_ptr<X1> const& __X1,
      boost::shared_ptr<RepType> __M,
      size_type rowstart = 0,
      bool init = false,
      bool do_threshold = false,
      typename X1::value_type threshold = type_traits<typename RepType::value_type>::epsilon() )
{
    return vf::detail::LinearForm<X1, RepType, RepType>( __X1, __M, rowstart, init, do_threshold, threshold );
}




/// \cond detail
template<typename Args>
struct compute_form1_return
{
#if 1
    typedef typename boost::remove_reference<typename parameter::binding<Args, tag::test>::type>::type::element_type test_type;
    //typedef typename boost::remove_reference<typename parameter::binding<Args, tag::vector>::type>::type::element_type vector_type;
    typedef typename Backend<double>::vector_type vector_type;
    typedef vf::detail::LinearForm<test_type,
            vector_type,
            vector_type> type;
#else
    typedef typename parameter::value_type<Args, tag::test>::type test_type;
    typedef typename parameter::value_type<Args, tag::vector>::type vector_type;


    typedef vf::detail::LinearForm<test_type,vector_type,vector_type> type;
#endif
};
/// \endcond
//boost::is_convertible<mpl::_,boost::shared_ptr<FunctionSpaceBase> >
BOOST_PARAMETER_FUNCTION(
    ( typename compute_form1_return<Args>::type ), // 1. return type
    form1,                                       // 2. name of the function template
    tag,                                        // 3. namespace of tag types
    ( required                                  // 4. one required parameter, and
      ( test,             *( boost::is_convertible<mpl::_,boost::shared_ptr<FunctionSpaceBase> > ) ) )
    ( optional                                  //    four optional parameters, with defaults
      //( in_out( vector ),   *( detail::is_vector_ptr<mpl::_> ), backend()->newVector( _test=test ) )
      ( backend,          *, Feel::backend() )
      ( in_out( vector ),   *, backend->newVector( test ) )
      ( init,             *( boost::is_integral<mpl::_> ), false )
      ( do_threshold,     *( boost::is_integral<mpl::_> ), bool( false ) )
      ( threshold,        *( boost::is_floating_point<mpl::_> ), type_traits<double>::epsilon() )
      ( rowstart,         *( boost::is_integral<mpl::_> ), 0 )
    )
)
{
    //Feel::detail::ignore_unused_variable_warning(boost_parameter_enabler_argument);
    Feel::detail::ignore_unused_variable_warning( args );
    //return form( test, *vector, init, false, 1e-16 );
    return form( test, vector, rowstart, init, do_threshold, threshold );
} // form

BOOST_PARAMETER_FUNCTION(
    ( typename compute_form1_return<Args>::type ), // 1. return type
    lform,                                       // 2. name of the function template
    tag,                                        // 3. namespace of tag types
    ( required                                  // 4. one required parameter, and
      ( test,             *( boost::is_convertible<mpl::_,boost::shared_ptr<FunctionSpaceBase> > ) )
      ( in_out( vector ),   *(Feel::detail::is_vector_ptr<mpl::_> ) )
        ) // required
    ( optional                                  //    four optional parameters, with defaults
      ( init,             *( boost::is_integral<mpl::_> ), false )
      ( do_threshold,     *( boost::is_integral<mpl::_> ), bool( false ) )
      ( threshold,        *( boost::is_floating_point<mpl::_> ), type_traits<double>::epsilon() )
      ( rowstart,         *( boost::is_integral<mpl::_> ), 0 )
    )
)
{
    //return form( test, *vector, init, false, 1e-16 );
    return form( test, *vector, rowstart, init, do_threshold, threshold );
} // form

/// \cond detail
template<typename Args, typename T>
struct compute_form2_return
{};

template<typename Args>
struct compute_form2_return<Args, mpl::false_>
{
    typedef typename parameter::value_type<Args, tag::test>::type::element_type::value_type value_type;
    typedef vf::detail::BilinearForm<typename parameter::value_type<Args, tag::test>::type::element_type,
            typename parameter::value_type<Args, tag::trial>::type::element_type,
            //typename parameter::value_type<Args, tag::matrix>::type::element_type,
            VectorUblas<value_type> > type;
};
template<typename Args>
struct compute_form2_return<Args, mpl::true_>
{
    typedef typename parameter::value_type<Args, tag::test>::type::element_type::value_type value_type;
    typedef vf::detail::BilinearForm<typename parameter::value_type<Args, tag::test>::type::element_type,
            typename parameter::value_type<Args, tag::test>::type::element_type,
            //typename parameter::value_type<Args, tag::vector>::type::element_type,
            VectorUblas<value_type> > type;
};
/// \endcond

BOOST_PARAMETER_FUNCTION( ( typename compute_form2_return<Args,mpl::bool_<boost::is_same<typename parameter::value_type<Args, tag::trial>::type, boost::parameter::void_>::value> >::type ), // 1. return type
                          form2,                                       // 2. name of the function template
                          tag,                                        // 3. namespace of tag types
                          ( required                                  // 4. one required parameter, and
                            ( test,             * )
                            ( trial,            * )
                          ) // required
                          (deduced
                           ( optional                                  //    four optional parameters, with defaults
                             ( init,             *( boost::is_integral<mpl::_> ), false )
                             ( pattern,          *( boost::is_integral<mpl::_> ), size_type( Pattern::COUPLED ) )
                             ( backend,          *, Feel::backend() )
                             ( in_out( matrix ),   *(boost::is_convertible<mpl::_, boost::shared_ptr<MatrixSparse<double>>>), backend->newMatrix( _test=test, _trial=trial, _pattern=pattern ) )
                             ( rowstart,         *( boost::is_integral<mpl::_> ), 0 )
                             ( colstart,         *( boost::is_integral<mpl::_> ), 0 )
                               ) // optional
                              ) // deduced
                        )
{
    Feel::detail::ignore_unused_variable_warning( args );
    //return form( test, trial, *matrix, init, false, 1e-16, pattern );
    //if (!matrix) matrix.reset( backend()->newMatrix( _trial=trial, _test=test ) );
    bool do_threshold = false;
    double threshold = 1e-16;
    return form( test, trial, matrix, rowstart, colstart, init, do_threshold, threshold, pattern );
    //return form( test, trial, *matrix, init, false, threshold, pattern );
    //return form( test, trial, *matrix, init, false, threshold, 0 );
} //


#if 0
BOOST_PARAMETER_FUNCTION(
    ( typename compute_form2_return<Args,mpl::bool_<boost::is_same<typename parameter::value_type<Args, tag::trial>::type, boost::parameter::void_>::value> >::type ), // 1. return type
    blform,                                       // 2. name of the function template
    tag,                                        // 3. namespace of tag types
    ( required                                  // 4. one required parameter, and
      ( test,             *( boost::is_convertible<mpl::_,boost::shared_ptr<FunctionSpaceBase> > ) )
      ( trial,            *( boost::is_convertible<mpl::_,boost::shared_ptr<FunctionSpaceBase> > ) )
      ( in_out( matrix ),   *(Feel::detail::is_matrix_ptr<mpl::_> ) ) ) // required
    ( optional                                  //    four optional parameters, with defaults
      ( init,             *( boost::is_integral<mpl::_> ), false )
      ( do_threshold,     *( boost::is_integral<mpl::_> ), bool( false ) )
      ( threshold,        *( boost::is_floating_point<mpl::_> ), type_traits<double>::epsilon() )
      ( pattern,          *( boost::is_integral<mpl::_> ), size_type( Pattern::COUPLED ) )
    )
)
{
    return form( test, trial, *matrix, init, do_threshold, threshold, pattern );
    //return form( test, trial, *matrix, init, false, 1e-16, pattern );
} //
#endif



} // Feel
#endif /* __Form_H */
