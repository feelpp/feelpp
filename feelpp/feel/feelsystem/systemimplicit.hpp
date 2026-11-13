/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 *
 * SPDX-FileCopyrightText: 2012-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
   \file systemimplicit.hpp
   \brief Legacy implicit system compatibility wrapper.
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-01-04
 */

#ifndef __SystemImplicit_H
#define __SystemImplicit_H 1

#include <feel/feelsystem/systemexplicit.hpp>

namespace Feel
{
/**
 * \class SystemImplicit
 * \brief describes an implicit system
 *
 * @author Christophe Prud'homme
 * @see Feel::SystemExplicit
 */
template<typename SpaceType>
class SystemImplicit : SystemExplicit<SpaceType>
{
    /** @brief Base legacy explicit system type. */
    typedef SystemExplicit<SpaceType> super;
public:


    /** @name Constants
     */
    //@{


    //@}

    /** @name Typedefs
     */
    //@{

    /** @brief Self type. */
    typedef SystemImplicit<SpaceType> system_type;

    /** @brief Scalar value type. */
    typedef typename super::value_type value_type;
    /** @brief Function-space type. */
    typedef typename super::functionspace_type functionspace_type;
    /** @brief Legacy function-space argument type. */
    typedef typename super::functionspace_type functionspace_ptrtype;
    /** @brief Finite-element field type. */
    typedef typename super::element_type element_type;

    /** @brief Backend type used by the legacy implicit system. */
    typedef Backend<value_type> backend_type;
    /** @brief Shared backend pointer type. */
    typedef std::shared_ptr<backend_type> backend_ptrtype;

    /*matrix*/
    /** @brief Sparse matrix type. */
    typedef typename backend_type::sparse_matrix_type sparse_matrix_type;
    /** @brief Sparse matrix pointer type. */
    typedef typename backend_type::sparse_matrix_ptrtype sparse_matrix_ptrtype;
    /* vector */
    /** @brief Vector type. */
    typedef typename backend_type::vector_type vector_type;
    /** @brief Vector pointer type. */
    typedef typename backend_type::vector_ptrtype vector_ptrtype;


    //@}

    /** @name Constructors, destructor
     */
    //@{

    /** @brief Construct from a function space and variables map. */
    SystemImplicit( functionspace_ptrtype const& Xh, po::variables_map const& vm );

    //! copy constructor
    SystemImplicit( SystemImplicit const & si );

    //! destructor
    ~SystemImplicit() {}

    //@}

    /** @name Operator overloads
     */
    //@{

    //! copy operator
    SystemImplicit& operator=( SystemImplicit const & o )
    {
        if ( this != &o )
        {
            super::operator=( o );
            M_backend = o.M_backend;
        }

        return *this;
    }
    //@}

    /** @name Accessors
     */
    //@{

    //! set the backend
    backend_ptrtype& backend()
    {
        return M_backend;
    }

    //! \return the backend
    backend_ptrtype const& backend() const
    {
        return M_backend;
    }

    //@}

    /** @name  Mutators
     */
    //@{

    //! set the backend
    void backend( backend_ptrtype const& b )
    {
        M_backend = b;
    }

    //@}

    /** @name  Methods
     */
    //@{


    //@}



protected:

    /** @brief Backend used by derived legacy systems. */
    backend_ptrtype M_backend;


private:

};
template<typename SpaceType>
SystemImplicit<SpaceType>::SystemImplicit( functionspace_ptrtype const& Xh,
        po::variables_map const& vm )
    :
    super( Xh, vm ),
    M_backend( backend_type::build( vm ) )
{}
template<typename SpaceType>
SystemImplicit<SpaceType>::SystemImplicit( SystemImplicit const& sil )
    :
    super( sil ),
    M_backend( sil.M_backend )
{}


} // Feel
#endif /* __SystemImplicit_H */
