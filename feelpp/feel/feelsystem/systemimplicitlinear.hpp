/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 *
 * SPDX-FileCopyrightText: 2012-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
   \file systemimplicitlinear.hpp
   \brief Legacy linear implicit system compatibility wrapper.
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-01-04
 */
#ifndef __SystemImplicitLinear_H
#define __SystemImplicitLinear_H 1

#include <feel/feelalg/backend.hpp>
#include <feel/feelsystem/systemimplicit.hpp>

namespace Feel
{
/**
 * \class SystemImplicitLinear
 * \brief Legacy implicit linear system with lhs/rhs storage.
 *
 * @author Christophe Prud'homme
 * @see Feel::SystemImplicit
 */
template<typename SpaceType>
class SystemImplicitLinear : public SystemImplicit<SpaceType>
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

    /** @brief Self type. */
    typedef SystemImplicitLinear<SpaceType> system_type;

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

    /** @brief Construct from a function space and variables map. */
    SystemImplicitLinear( functionspace_ptrtype const& Xh, po::variables_map const& vm );
    /** @brief Copy constructor. */
    SystemImplicitLinear( SystemImplicitLinear const & sil );
    /** @brief Destructor. */
    ~SystemImplicitLinear() {}

    //@}

    /** @name Operator overloads
     */
    //@{

    /** @brief Copy assignment. */
    SystemImplicitLinear& operator=( SystemImplicitLinear const & o )
    {
        if ( this != &o )
        {
            super::operator=( o );


            M_lhs = o.M_lhs;
            M_rhs = o.M_rhs;
        }

        return *this;
    }
    //@}

    /** @name Accessors
     */
    //@{

    //! \return left hand side (lhs) matrix
    sparse_matrix_ptrtype const& lhs() const
    {
        return M_lhs;
    }

    //! \return left hand side (lhs) matrix
    sparse_matrix_ptrtype& lhs()
    {
        return M_lhs;
    }

    //! \return right hand side (rhs) vector
    vector_ptrtype const& rhs() const
    {
        return M_rhs;
    }

    //! \return right hand side (rhs) vector
    vector_ptrtype& rhs()
    {
        return M_rhs;
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
     * solve the linear system
     */
    void solve( element_type& u )
    {
        vector_ptrtype U( this->backend->newVector( _test=this->functionSpace() ) );
        this->backend->solve( M_lhs, M_lhs, U, M_rhs );
        u = *U;
    }

    //@}



protected:

    /** @brief Legacy lhs matrix storage. */
    sparse_matrix_ptrtype M_lhs;
    /** @brief Legacy rhs vector storage. */
    vector_ptrtype M_rhs;

private:



};
template<typename SpaceType>
SystemImplicitLinear<SpaceType>::SystemImplicitLinear( functionspace_ptrtype const& Xh,
        po::variables_map const& vm )
    :
    super( Xh, vm ),
    M_lhs( this->backend()->newMatrix( Xh, Xh ) ),
    M_rhs( this->backend()->newVector( Xh ) )
{}
template<typename SpaceType>
SystemImplicitLinear<SpaceType>::SystemImplicitLinear( SystemImplicitLinear const& sil )
    :
    super( sil ),
    M_lhs( sil.M_lhs ),
    M_rhs( sil.M_rhs )
{}


}
#endif /* __SystemImplicitLinear_H */
