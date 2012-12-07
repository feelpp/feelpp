/* -*- mode: c++ -*-

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2010-03-03

  Copyright (C) 2010 Université Joseph Fourier (Grenoble I)

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
   \file operatorlinearparametrized.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2010-03-03
 */
#ifndef __OperatorLinearParametrized_H
#define __OperatorLinearParametrized_H 1

#include <feel/feelcrb/parameterspace.hpp>

namespace Feel
{
/**
 * \class OperatorLinearParametrized
 * \brief An interface for linear parametrized operators
 *
 * @author Christophe Prud'homme
 * @see
 */
template<class DomainSpace, class DualImageSpace>
class OperatorLinearParametrized : public OperatorLinear<DomainSpace, DualImageSpace>
{
    typedef OperatorLinearParametrized<DomainSpace,DualImageSpace> super;
public:


    /** @name Constants
     */
    //@{


    //@}

    /** @name Typedefs
     */
    //@{
    // -- TYPEDEFS --
    typedef OperatorLinearParametrized<DomainSpace, DualImageSpace> this_type;
    typedef OperatorLinear<DomainSpace, DualImageSpace> super_type;

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

    /* parameter space */
    typedef ParameterSpace<ParameterSpaceDimension> parameterspace_type;
    typedef boost::shared_ptr<parameterspace_type> parameterspace_ptrtype;
    typedef parameterspace_type::element_type parameter_type;
    typedef parameterspace_type::element_ptrtype parameter_ptrtype;
    typedef parameterspace_type::sampling_type sampling_type;
    typedef parameterspace_type::sampling_ptrtype sampling_ptrtype;

    typedef Eigen::VectorXd theta_vector_type;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    OperatorLinearParametrized()
        :
        super_type()
    {}

    //! copy constructor
    OperatorLinearParametrized( OperatorLinearParametrized const & olp, bool deep_copy = false )
        :
        super_type( olp, deep_copy )
    {}

    /**
     * Constructor from domain and image space
     * \param domainSpace
     * \param dualImageSpace
     * \param backend associated linear algebra backend
     */
    OperatorLinearParametrized( domain_space_ptrtype     domainSpace,
                                dual_image_space_ptrtype dualImageSpace,
                                backend_ptrtype          backend )
        :
        super_type( domainSpace, dualImageSpace, backend )
    {
    }

    //! destructor
    ~OperatorLinearParametrized()
    {}

    void
    init( domain_space_ptrtype     domainSpace,
          dual_image_space_ptrtype dualImageSpace,
          backend_ptrtype          backend )
    {
        super::init( domainSpace, dualImageSpace, backend );
    }
    //@}

    /** @name Operator overloads
     */
    //@{

    //! copy operator
    OperatorLinearParametrized& operator=( OperatorLinearParametrized const & o )
    {
        if ( this != &o )
        {
        }

        return *this;
    }
    //@}

    /** @name Accessors
      */
    //@{

    //! return the parameter space
    parameterspace_ptrtype parameterSpace() const
    {
        return M_Dmu;
    }

    /**
     * \brief return the coefficient vector
     */
    theta_vector_type const& thetaq() const
    {
        return M_thetaq;
    }

    /**
     * \brief return the coefficient vector \p q component
     *
     */
    value_type thetaq( int q ) const
    {
        return M_thetaq( q );
    }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    /**
     * \brief compute the theta coefficient for both bilinear and linear form
     * \param mu parameter to evaluate the coefficients
     */
    virtual theta_vector_type computeThetaq( parameter_type const& mu ) = 0;

    // fill underlying matrix
    template<class ExprT>
    this_type& add( int q, ExprT const& e )
    {
        //         M_matrix->clear();
        form2( this->domainSpace(),
               this->dualImageSpace(),
               M_Aq[q],
               _init=true

             ) = e;
        return *this;
    }

    // add to underlying matrix
    template<class ExprT>
    this_type& operator+=( ExprT const& e )
    {
        form2( this->domainSpace(),
               this->dualImageSpace(),
               M_matrix,
               _init=false ) += e;
        return *this;
    }

    /**
     * \brief update the model wrt \p mu
     */
    offline_merge_type update( parameter_type const& mu )
    {
        this->computeThetaq( mu );
        return offlineMerge( mu );
    }

    /**
     * Sum the affine decomposition for the parameter mu, need to call
     * update( mu ) before calling the merge( mu )
     */
    sparse_matrix_type merge( parameter_type const& mu );

    //@}



protected:

    /**
     * matrix storing the parameter independent matrices (affine
     * decomposition)
     */
    std::vector<sparse_matrix_ptrtype> M_Aq;

    /**
     * parameter space
     */
    parameterspace_ptrtype M_Dmu;

    //! coefficients of the Aq matrices
    theta_vector_type M_thetaq;
};

template<class DomainSpace, class DualImageSpace>
typename OperatorLinearParametrized<DomainSpace,DualImageSpace>::sparse_matrix_ptrtype
OperatorLinearParametrized<DomainSpace,DualImageSpace>::merge( parameter_type const& mu );
{
    sparse_matrix_ptrtype A( M_backend->newMatrix( domainSpace(), dualImageSpace() ) )

    A->close();
    *A = *M_Aq[0];
    A->scale( this->thetaq( 0 ) );

    for ( int q = 1; q < Qa(); ++q )
    {
        A->addMatrix( this->thetaq( q ), M_Aq[q] );
    }

    return A;
}
} // Feel
#endif /* __OperatorLinearParametrized_H */
