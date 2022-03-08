/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2013-04-26

  Copyright (C) 2013-2016 Feel++ Consortium

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
   \file operatorlinearfree.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \author Stephane Veys <stephane.veys@imag.fr>
   \date 2013-04-26
 */
#ifndef FEELPP_DISCR_OPERATORLINEARFREE_H
#define FEELPP_DISCR_OPERATORLINEARFREE_H 1

#include <boost/ref.hpp>
#include <boost/next_prior.hpp>
#include <boost/type_traits.hpp>
#include <boost/tuple/tuple.hpp>
#if BOOST_VERSION >= 104700
#include <boost/math/special_functions/nonfinite_num_facets.hpp>
#endif
#include <feel/feeldiscr/operatorlinear.hpp>
#include <feel/feelvf/vf.hpp>
namespace Feel
{

/**
 * \class OperatorLinearFree
 * \brief Linear Operator Free between function spaces, represented by a matrix
 * but this matrix is not automatically assembled
 */
template<class DomainSpace, class DualImageSpace, typename ExprType>
class OperatorLinearFree : public OperatorLinear<DomainSpace, DualImageSpace>
{
    typedef OperatorLinear<DomainSpace,DualImageSpace> super;
public:

    typedef ExprType expr_type;

    typedef typename super::domain_space_type domain_space_type;
    typedef typename super::dual_image_space_type  dual_image_space_type;
    typedef typename super::domain_space_ptrtype domain_space_ptrtype;
    typedef typename super::dual_image_space_ptrtype  dual_image_space_ptrtype;

    typedef typename super::backend_type backend_type;
    typedef typename super::backend_ptrtype backend_ptrtype;

    typedef typename super::matrix_type matrix_type;
    typedef std::shared_ptr<matrix_type> matrix_ptrtype;

    typedef FsFunctionalLinear<DualImageSpace> image_element_type;
    typedef typename super::domain_element_type domain_element_type;
    typedef typename super::dual_image_element_type dual_image_element_type;

    typedef typename super::domain_element_range_type domain_element_range_type;
    typedef typename super::domain_element_slice_type domain_element_slice_type;
    typedef typename super::dual_image_element_range_type dual_image_element_range_type;
    typedef typename super::dual_image_element_slice_type dual_image_element_slice_type;

    typedef OperatorLinearFree<DomainSpace, DualImageSpace, ExprType> this_type;
    typedef std::shared_ptr<this_type> this_ptrtype;

    OperatorLinearFree()
        :
        super()
    {}

    OperatorLinearFree (domain_space_ptrtype     domainSpace,
                        dual_image_space_ptrtype dualImageSpace,
                        backend_ptrtype          backend,
                        expr_type                expr,
                        size_type                pattern=Pattern::COUPLED )
        :
        super( domainSpace , dualImageSpace , backend , false , pattern ),
        M_backend( backend ),
        M_expr( expr ),
        M_pattern( pattern )
    {}

    ~OperatorLinearFree() override {}


    void
    init( domain_space_ptrtype     domainSpace,
          dual_image_space_ptrtype dualImageSpace,
          backend_ptrtype          backend,
          bool buildMatrix = false ,
          size_type pattern = Pattern::COUPLED ) override
    {
        this->setDomainSpace( domainSpace );
        this->setDualImageSpace( dualImageSpace );
        M_backend = backend;
        M_pattern = pattern;
    }


    expr_type expr()
    {
        return M_expr;
    }

    void
    apply( const domain_element_type& de,
           image_element_type&        ie ) const override
    {

        auto matrix = M_backend->newMatrix( _test=this->dualImageSpace() ,
                                            _trial=this->domainSpace(),
                                            _pattern=M_pattern );


        form2( _trial=this->domainSpace(),
               _test=this->dualImageSpace(),
               _matrix=matrix,
               _pattern=M_pattern) = M_expr;

        matrix->close();

        vector_ptrtype _v1( M_backend->newVector( _test=de.functionSpace() ) );
        *_v1 = de;_v1->close();
        vector_ptrtype _v2( M_backend->newVector( _test=ie.space() ) );
        M_backend->prod( matrix, _v1, _v2 );
        ie.container() = *_v2;
    }


    double
    energy( const typename domain_space_type::element_type & de,
            const typename dual_image_space_type::element_type & ie ) const override
    {

        auto matrix = M_backend->newMatrix( _test=this->dualImageSpace() ,
                                            _trial=this->domainSpace(),
                                            _pattern=M_pattern );

        form2( _trial=this->domainSpace(),
               _test=this->dualImageSpace(),
               _matrix=matrix,
               _pattern=M_pattern) = M_expr;

        matrix->close();

        vector_ptrtype _v1( M_backend->newVector( _test=de.functionSpace() ) );
        *_v1 = de;_v1->close();
        vector_ptrtype _v2( M_backend->newVector( _test=ie.functionSpace() ) );
        *_v2 = ie;
        vector_ptrtype _v3( M_backend->newVector( _test=ie.functionSpace() ) );
        M_backend->prod( matrix, _v1, _v3 );
        return inner_product( _v2, _v3 );
    }


    void
    apply( const typename domain_space_type::element_type & de,
           typename dual_image_space_type::element_type & ie ) override
    {
        auto matrix = M_backend->newMatrix( _test=this->dualImageSpace() ,
                                            _trial=this->domainSpace(),
                                            _pattern=M_pattern );

        form2( _trial=this->domainSpace(),
               _test=this->dualImageSpace(),
               _matrix=matrix,
               _pattern=M_pattern) = M_expr;

        matrix->close();

        vector_ptrtype _v1( M_backend->newVector( _test=de.functionSpace() ) );
        *_v1 = de;_v1->close();
        vector_ptrtype _v2( M_backend->newVector( _test=ie.functionSpace() ) );
        M_backend->prod( matrix, _v1, _v2 );
        ie.container() = *_v2;
    }


    void
    apply( const domain_element_range_type & de,
           typename dual_image_space_type::element_type & ie ) override
    {
        auto matrix = M_backend->newMatrix( _test=this->dualImageSpace() ,
                                            _trial=this->domainSpace(),
                                            _pattern=M_pattern );

        form2( _trial=this->domainSpace(),
               _test=this->dualImageSpace(),
               _matrix=matrix,
               _pattern=M_pattern) = M_expr;

        matrix->close();

        vector_ptrtype _v1( M_backend->newVector( _test=de.functionSpace() ) );
        *_v1 = de;_v1->close();
        vector_ptrtype _v2( M_backend->newVector( _test=ie.functionSpace() ) );
        M_backend->prod( matrix, _v1, _v2 );
        ie.container() = *_v2;
    }

    void
    apply( const typename domain_space_type::element_type & de,
           dual_image_element_range_type & ie ) override
    {

        auto matrix = M_backend->newMatrix( _test=this->dualImageSpace() ,
                                            _trial=this->domainSpace(),
                                            _pattern=M_pattern );

        form2( _trial=this->domainSpace(),
               _test=this->dualImageSpace(),
               _matrix=matrix,
               _pattern=M_pattern) = M_expr;

        matrix->close();

        vector_ptrtype _v1( M_backend->newVector( _test=de.functionSpace() ) );
        *_v1 = de;_v1->close();
        vector_ptrtype _v2( M_backend->newVector( _test=ie.functionSpace() ) );
        M_backend->prod( matrix, _v1, _v2 );
        ie.container() = *_v2;
    }

    void
    apply( const domain_element_range_type & de,
           dual_image_element_range_type & ie ) override
    {

        auto matrix = M_backend->newMatrix( _test=this->dualImageSpace() ,
                                            _trial=this->domainSpace(),
                                            _pattern=M_pattern );

        form2( _trial=this->domainSpace(),
               _test=this->dualImageSpace(),
               _matrix=matrix,
               _pattern=M_pattern) = M_expr;

        matrix->close();

        vector_ptrtype _v1( M_backend->newVector( _test=de.functionSpace() ) );
        *_v1 = de;_v1->close();
        vector_ptrtype _v2( M_backend->newVector( _test=ie.functionSpace() ) );
        M_backend->prod( matrix, _v1, _v2 );
        ie.container() = *_v2;
    }

    void
    apply( const domain_element_slice_type & de,
           typename dual_image_space_type::element_type & ie ) override
    {
        auto matrix = M_backend->newMatrix( _test=this->dualImageSpace() ,
                                            _trial=this->domainSpace(),
                                            _pattern=M_pattern );

        form2( _trial=this->domainSpace(),
               _test=this->dualImageSpace(),
               _matrix=matrix,
               _pattern=M_pattern) = M_expr;

        matrix->close();

        vector_ptrtype _v1( M_backend->newVector( _test=de.functionSpace() ) );
        *_v1 = de;_v1->close();
        vector_ptrtype _v2( M_backend->newVector( _test=ie.functionSpace() ) );
        M_backend->prod( matrix, _v1, _v2 );
        ie.container() = *_v2;
    }


    void
    apply( const typename domain_space_type::element_type & de,
           dual_image_element_slice_type & ie ) override
    {
        auto matrix = M_backend->newMatrix( _test=this->dualImageSpace() ,
                                            _trial=this->domainSpace(),
                                            _pattern=M_pattern );

        form2( _trial=this->domainSpace(),
               _test=this->dualImageSpace(),
               _matrix=matrix,
               _pattern=M_pattern) = M_expr;

        matrix->close();

        vector_ptrtype _v1( M_backend->newVector( _test=de.functionSpace() ) );
        *_v1 = de;_v1->close();
        vector_ptrtype _v2( M_backend->newVector( _test=ie.functionSpace() ) );
        M_backend->prod( matrix, _v1, _v2 );
        ie.container() = *_v2;
    }

    void
    apply( /*const*/ domain_element_slice_type /*&*/ de,
                     dual_image_element_slice_type /*&*/ ie ) override
    {
        auto matrix = M_backend->newMatrix( _test=this->dualImageSpace() ,
                                            _trial=this->domainSpace(),
                                            _pattern=M_pattern );

        form2( _trial=this->domainSpace(),
               _test=this->dualImageSpace(),
               _matrix=matrix,
               _pattern=M_pattern) = M_expr;

        matrix->close();

        vector_ptrtype _v1( M_backend->newVector( _test=de.functionSpace() ) );
        *_v1 = de;_v1->close();
        vector_ptrtype _v2( M_backend->newVector( _test=ie.functionSpace() ) );
        M_backend->prod( matrix, _v1, _v2 );
        ie.container() = *_v2;
    }



    void
    apply( const domain_element_range_type & de,
           dual_image_element_slice_type & ie ) override
    {
        auto matrix = M_backend->newMatrix( _test=this->dualImageSpace() ,
                                            _trial=this->domainSpace(),
                                            _pattern=M_pattern );

        form2( _trial=this->domainSpace(),
               _test=this->dualImageSpace(),
               _matrix=matrix,
               _pattern=M_pattern) = M_expr;

        matrix->close();

        vector_ptrtype _v1( M_backend->newVector( _test=de.functionSpace() ) );
        *_v1 = de;_v1->close();
        vector_ptrtype _v2( M_backend->newVector( _test=ie.functionSpace() ) );
        M_backend->prod( matrix, _v1, _v2 );
        ie.container() = *_v2;
    }

    void
    apply( const domain_element_slice_type & de,
           dual_image_element_range_type & ie ) override
    {
        auto matrix = M_backend->newMatrix( _test=this->dualImageSpace() ,
                                            _trial=this->domainSpace(),
                                            _pattern=M_pattern );

        form2( _trial=this->domainSpace(),
               _test=this->dualImageSpace(),
               _matrix=matrix,
               _pattern=M_pattern) = M_expr;

        matrix->close();

        vector_ptrtype _v1( M_backend->newVector( _test=de.functionSpace() ) );
        *_v1 = de;_v1->close();
        vector_ptrtype _v2( M_backend->newVector( _test=ie.functionSpace() ) );
        M_backend->prod( matrix, _v1, _v2 );
        ie.container() = *_v2;
    }

    //! apply the inverse of the operator: \f$de = O^{-1} ie\f$
    void
    applyInverse( domain_element_type&      de,
                  const image_element_type& ie ) override
    {
        auto matrix = M_backend->newMatrix( _test=this->dualImageSpace() ,
                                            _trial=this->domainSpace(),
                                            _pattern=M_pattern );

        form2( _trial=this->domainSpace(),
               _test=this->dualImageSpace(),
               _matrix=matrix,
               _pattern=M_pattern) = M_expr;

        matrix->close();

        vector_ptrtype _v1( M_backend->newVector( _test=de.functionSpace() ) );
        vector_ptrtype _v2( M_backend->newVector( _test=ie.space() ) );
        *_v2 = ie.container();
        M_backend->solve( matrix, matrix, _v1, _v2 );
        de = *_v1;

    }


    this_type& operator=( this_type const& m )
    {
        M_backend = m.M_backend;
        M_expr = m.M_expr;
        M_pattern = m.M_pattern;
        return *this;
    }


    // retrieve underlying matrix
    void matPtr( matrix_ptrtype & matrix ) override
    {
        form2( _trial=this->domainSpace(),
               _test=this->dualImageSpace(),
               _matrix=matrix,
               _pattern=M_pattern) = M_expr;

        matrix->close();
    }

private:
    backend_ptrtype M_backend;
    expr_type M_expr;
    size_type M_pattern;
};//OperatorLinearFree


template <typename ... Ts>
auto opLinearFree( Ts && ... v )
{
    auto args = NA::make_arguments( std::forward<Ts>(v)... );
    auto && domainSpace = args.get(_domainSpace);
    auto && imageSpace = args.get(_imageSpace);
    auto && expr = args.get(_expr);
    auto && backend = args.get_else_invocable(_backend, [](){ return Feel::backend(); } );
    size_type pattern = args.get_else( _pattern, Pattern::COUPLED );

    using domain_space_type = Feel::remove_shared_ptr_type<std::remove_pointer_t<std::decay_t<decltype(domainSpace)>>>;
    using image_space_type = Feel::remove_shared_ptr_type<std::remove_pointer_t<std::decay_t<decltype(imageSpace)>>>;
    using expr_type = std::decay_t<decltype(expr)>;
    return std::make_shared<OperatorLinearFree<domain_space_type, image_space_type, expr_type>>( domainSpace,imageSpace,backend,expr );
}


}//namespace Feel
#endif
