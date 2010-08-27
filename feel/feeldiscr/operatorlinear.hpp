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
   \file operatorlinear.hpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2006-11-16
 */

#ifndef _OPERATORLINEAR_HPP_
#define _OPERATORLINEAR_HPP_

#include <feel/feeldiscr/operator.hpp>
#include <feel/feelvf/vf.hpp>

namespace Feel
{
/**
 * \class OperatorLinear
 * \brief Linear Operator between function spaces, represented by a matrix
 */
template<class DomainSpace, class DualImageSpace>
class OperatorLinear : public Operator<DomainSpace, DualImageSpace>
{
    typedef OperatorLinear<DomainSpace,DualImageSpace> super;
public:

    // -- TYPEDEFS --
    typedef OperatorLinear<DomainSpace, DualImageSpace> this_type;
    typedef Operator<DomainSpace, DualImageSpace> super_type;

    typedef typename super::domain_space_type domain_space_type;
    typedef typename super::dual_image_space_type  dual_image_space_type;
    typedef typename super::domain_space_ptrtype domain_space_ptrtype;
    typedef typename super::dual_image_space_ptrtype  dual_image_space_ptrtype;
    typedef typename domain_space_type::element_type domain_element_type;

    typedef typename super::backend_type backend_type;
    typedef typename super::backend_ptrtype backend_ptrtype;
    typedef typename backend_type::sparse_matrix_type matrix_type;
    typedef typename backend_type::vector_type vector_type;
    typedef typename backend_type::vector_ptrtype vector_ptrtype;
    typedef boost::shared_ptr<matrix_type> matrix_ptrtype;

    template<typename T, typename Storage>
    struct domain_element: public super::domain_space_type::template Element<T,Storage> {};

    typedef FsFunctionalLinear<DualImageSpace> image_element_type;

    OperatorLinear()
        :
        super_type(),
#if defined( HAVE_PETSC_H )
        M_backend( backend_type::build( BACKEND_PETSC ) ),
#else
        M_backend( backend_type::build( BACKEND_GMM ) ),
#endif
        M_matrix()
    {}
    OperatorLinear( OperatorLinear const& ol, bool deep_copy = false )
        :
        super_type( ol ),
        M_backend( ol.M_backend ),
        M_matrix()
    {
        if ( deep_copy )
            {
                M_matrix = new matrix_type( *ol.M_matrix );
            }
        else
            {
                M_matrix = ol.M_matrix;
            }
    }
    OperatorLinear( domain_space_ptrtype     domainSpace,
                    dual_image_space_ptrtype dualImageSpace,
                    backend_ptrtype          backend ) :
        super_type( domainSpace, dualImageSpace ),
        M_backend( backend ),
        M_matrix( M_backend->newMatrix( domainSpace, dualImageSpace ) )
    {
    }

    ~OperatorLinear() {}

    void
    init( domain_space_ptrtype     domainSpace,
          dual_image_space_ptrtype dualImageSpace,
          backend_ptrtype          backend )
    {
        this->setDomainSpace( domainSpace );
        this->setDualImageSpace( dualImageSpace );
        M_backend = backend;
        M_matrix = M_backend->newMatrix( domainSpace, dualImageSpace );


    }

    // apply the operator: ie := Op de
    template<typename Storage>
    void
    apply( const domain_element<typename domain_element_type::value_type, Storage>& de,
           image_element_type&        ie ) const
    {
        if ( ! M_matrix->closed() )
        {
            M_matrix->close();
        }
        vector_ptrtype _v1( M_backend->newVector( de.map() ) );
        *_v1 = de;
        vector_ptrtype _v2( M_backend->newVector( ie.space()->map() ) );
        M_backend->prod( M_matrix, _v1, _v2 );
        ie.container() = *_v2;
    }

    void
    apply( const domain_element_type& de,
           image_element_type&        ie ) const
    {
        if ( ! M_matrix->closed() )
        {
            M_matrix->close();
        }
        vector_ptrtype _v1( M_backend->newVector( de.map() ) );
        *_v1 = de;
        vector_ptrtype _v2( M_backend->newVector( ie.space()->map() ) );
        M_backend->prod( M_matrix, _v1, _v2 );
        ie.container() = *_v2;
    }


    //template <typename T1 = typename domain_space_type::element_type,
    //          typename T2 = typename dual_image_space_type::element_type >
    //void
    //applyBis( T1 & de, T2 & ie ) //const
    void
    apply( typename domain_space_type::element_type & de,
           typename dual_image_space_type::element_type & ie)
    {
        if ( ! M_matrix->closed() )
        {
            M_matrix->close();
        }
        vector_ptrtype _v1( M_backend->newVector( de.map() ) );
        *_v1 = de;
        //vector_ptrtype _v2( M_backend->newVector( ie.space()->map() ) );
        vector_ptrtype _v2( M_backend->newVector( ie.map() ) );
        M_backend->prod( M_matrix, _v1, _v2 );
        ie.container() = *_v2;
    }

    template <typename T1 = typename domain_space_type::element_type,
              typename T2 = typename dual_image_space_type::element_type >
    T2
    operator()(T1 & de)
    {
        T2 elt_image(this->dualImageSpace(),"oio");
        this->apply(de,elt_image);

        return elt_image;
    }


    //! apply the inverse of the operator: \f$de = O^{-1} ie\f$
    virtual void
    applyInverse( domain_element_type&      de,
                  const image_element_type& ie )
    {
        if ( ! M_matrix->closed() )
        {
            M_matrix->close();
        }
        vector_ptrtype _v1( M_backend->newVector( de.map() ) );
        vector_ptrtype _v2( M_backend->newVector( ie.space()->map() ) );
        *_v2 = ie.container();
        M_backend->solve( M_matrix, M_matrix, _v1, _v2 );
        de = *_v1;

#if 0
        if ( !M_backend->converged() )
        {
            std::cerr << "[OperatorLinear::applyInverse] "
                      << "solver failed to converge" << std::endl;
        }
#endif
    }

    // fill underlying matrix
    template<class ExprT>
    this_type& operator=( ExprT const& e )
    {
        //         M_matrix->clear();
        form2( _trial=this->domainSpace(),
               _test=this->dualImageSpace(),
               _matrix=M_matrix,
               _init=true
              ) = e;
        return *this;
    }

    this_type& operator=( this_type const& m )
    {
        M_backend = m.M_backend;
        M_matrix->zero();
        M_matrix->addMatrix( 1.0, m.M_matrix );

        return *this;
    }

    // add to underlying matrix
    template<class ExprT>
    this_type& operator+=( ExprT const& e )
    {
        form2( _trial=this->domainSpace(),
               _test=this->dualImageSpace(),
               _matrix=M_matrix ) += e;
        return *this;
    }

    void close()
    {
        if ( ! M_matrix->closed() )
            {
                M_matrix->close();
            }

    }
    // retrieve underlying matrix
    matrix_type& mat()
    {
        return *M_matrix;
    }

    // retrieve underlying matrix
    matrix_type const& mat() const
    {
        return *M_matrix;
    }

    // retrieve underlying matrix
    matrix_ptrtype const& matPtr() const
    {
        return M_matrix;
    }

    // retrieve underlying matrix
    matrix_ptrtype& matPtr()
    {
        return M_matrix;
    }

    template<typename T>
    OperatorLinear& add( T const& scalar, OperatorLinear const& ol )
    {
        this->close();
        M_matrix->addMatrix( scalar, *ol.M_matrix );
        return *this;
    }

    backend_ptrtype& backend()
    {
        return M_backend;
    }

    template<typename T>
    OperatorLinear& add( T const& scalar, boost::shared_ptr<OperatorLinear> ol )
    {
        this->close();
        this->add(scalar, *ol);
        return *this;
    }

    OperatorLinear& operator+( boost::shared_ptr<OperatorLinear> ol )
    {
        this->close();
        this->add(1.0, *ol);
        return *this;
    }
    OperatorLinear& operator+( OperatorLinear const& ol )
    {
        this->close();
        this->add(1.0, ol);
        return *this;
    }

private:

    backend_ptrtype M_backend;
    matrix_ptrtype M_matrix;


}; // class Operator

} // Feel

#endif /* _OPERATORLINEAR_HPP_ */
