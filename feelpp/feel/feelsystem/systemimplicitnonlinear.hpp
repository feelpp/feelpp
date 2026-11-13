/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 *
 * SPDX-FileCopyrightText: 2012-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
   \file systemimplicitnonlinear.hpp
   \brief Legacy nonlinear implicit system compatibility wrapper.
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-01-04
 */
#ifndef __SystemImplicitNonLinear_H
#define __SystemImplicitNonLinear_H 1

#include <feel/feelsystem/systemimplicit.hpp>


namespace Feel
{
/**
 * \class SystemImplicitNonLinear
 * \brief Describes nonlinear implicit systems
 *
 * @author Christophe Prud'homme
 * @see Feel::SystemImplicit
 */
template<typename SpaceType>
class SystemImplicitNonLinear : public SystemImplicit<SpaceType>
{
    /** @brief Base legacy implicit system type. */
    typedef SystemImplicit<SpaceType> super;
public:


    /** @name Constants
     */
    //@{


    //@}

    /** @name Typedefs
     */
    //@{

    /** @brief Self type used for callback binding. */
    typedef SystemImplicitNonLinear<SpaceType> self_type;
    /** @brief Self type. */
    typedef SystemImplicitNonLinear<SpaceType> system_type;

    /** @brief Scalar value type. */
    typedef typename super::value_type value_type;
    /** @brief Function-space type. */
    typedef typename super::functionspace_type functionspace_type;
    /** @brief Legacy function-space argument type. */
    typedef typename super::functionspace_type functionspace_ptrtype;
    /** @brief Finite-element field type. */
    typedef typename super::element_type element_type;

    /** @brief Backend type. */
    typedef typename super::backend_type backend_type;
    /** @brief Shared backend pointer type. */
    typedef typename super::backend_ptrtype backend_ptrtype;
    /** @brief Sparse matrix type. */
    typedef typename super::sparse_matrix_type sparse_matrix_type;
    /** @brief Sparse matrix pointer type. */
    typedef typename super::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    /** @brief Vector type. */
    typedef typename super::vector_type vector_type;
    /** @brief Vector pointer type. */
    typedef typename super::vector_ptrtype vector_ptrtype;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    //! default constructor
    SystemImplicitNonLinear( functionspace_ptrtype const& Xh, po::variables_map const& vm );
    //! copy constructor
    SystemImplicitNonLinear( SystemImplicitNonLinear const & );
    //! destructor
    ~SystemImplicitNonLinear();

    //@}

    /** @name Operator overloads
     */
    //@{

    //! copy operator
    SystemImplicitNonLinear& operator=( SystemImplicitNonLinear const & o )
    {
        if ( this != &o )
        {
            super::operator=( o );
            M_J = o.M_J;
            M_R = o.M_R;
        }

        return *this;
    }
    //@}

    /** @name Accessors
     */
    //@{

    //! \return the jacobian matrix
    sparse_matrix_ptrtype const& jacobian() const
    {
        return M_J;
    }

    //! \return the jacobian matrix
    sparse_matrix_ptrtype& jacobian()
    {
        return M_J;
    }

    //! \return the residual
    vector_ptrtype const& residual() const
    {
        return M_R;
    }

    //! \return the residual
    vector_ptrtype & residual()
    {
        return M_R;
    }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{

    //! function that update the residual at each non linear iteration
    virtual void updateResidual( const vector_ptrtype& X, vector_ptrtype& R ) = 0;

    //! function that update the jacobian at each non linear iteration
    virtual void updateJacobian( const vector_ptrtype& X, sparse_matrix_ptrtype& J ) = 0;

    //! default implementation that solves the nonlinear system
    void solve( element_type& u )
    {
        vector_ptrtype U( M_backend->newVector( this->functionSpace() ) );
        *U = u;

        this->updateResidual( U, M_R );
        this->updateJacobian( U, M_J );

        M_backend->nlSolve( M_J, U, M_R, 1e-10, 10 );
        u = *U;
    }
    //@}



protected:

    /** @brief Legacy Jacobian matrix storage. */
    sparse_matrix_ptrtype M_J;
    /** @brief Legacy residual vector storage. */
    vector_ptrtype M_R;

private:

};
template<typename SpaceType>
SystemImplicitNonLinear<SpaceType>::SystemImplicitNonLinear( functionspace_ptrtype const& Xh,
        po::variables_map const& vm )
    :
    super( Xh, vm ),
    M_J( M_backend->newMatrix( Xh, Xh ) ),
    M_R( M_backend->newVector( Xh ) )

{
    using namespace std::placeholders;
    M_backend->nlSolver()->residual = std::bind( &self_type::updateResidual, std::ref( *this ), _1, _2 );
    M_backend->nlSolver()->jacobian = std::bind( &self_type::updateJacobian, std::ref( *this ), _1, _2 );

}
template<typename SpaceType>
SystemImplicitNonLinear<SpaceType>::SystemImplicitNonLinear( SystemImplicitNonLinear const& sil )
    :
    super( sil ),
    M_J( sil.M_J ),
    M_R( sil.M_R )

{}


}
#endif /* __SystemImplicitNonLinear_H */
