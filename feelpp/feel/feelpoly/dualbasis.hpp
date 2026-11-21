// SPDX-FileCopyrightText: 2005-2006 EPFL
// SPDX-FileCopyrightText: 2026 University of Strasbourg
// SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
   \file dualbasis.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2005-10-10
 */
#ifndef __DualBasis_H
#define __DualBasis_H 1

#include <feel/feelpoly/order.hpp>

namespace Feel
{
/**
 * \class DualBasis
 * \brief basis of a space P' dual of some space P
 *
 * The element of this basis are functionals evaluated at a node set
 * It contains also the information about the dof table
 *
 * \ingroup Polynomial
 * @author Christophe Prud'homme
 */
template<typename Primal>
class DualBasis
{
public:


    /** @name Typedefs
     */
    //@{
    typedef DualBasis<Primal> self_type;

    typedef Primal primal_space_type;
    typedef typename primal_space_type::basis_type basis_type;

    //@}
    /** @name Constants
    */
    //@{
    static inline const uint16_type nDim = primal_space_type::nDim;
    /** Authoritative primal expansion order specification. */
    static constexpr int order_spec = orderSpecification<Primal>();
    /** Whether the primal carries runtime order. */
    static constexpr bool is_order_dynamic = order_spec == Dynamic;
    /** Whether the primal degree is fixed in its type. */
    static constexpr bool is_order_static = !is_order_dynamic;
    /** Legacy static geometry extent; runtime construction uses explicit order. */
    static constexpr uint16_type nOrder = is_order_dynamic ? 1 : order_spec;
    //@}
    /** @name Constructors, destructor
     */
    //@{


    DualBasis( primal_space_type const& primal )
        :
        M_primal( primal )
    {}
    DualBasis( DualBasis const & b )
        :
        M_primal( b.M_primal )
    {}
    ~DualBasis()
    {}

    //@}

    /** @name Operator overloads
     */
    //@{

    self_type const& operator=( self_type const& dual )
    {
        if ( this != &dual )
        {
            M_primal = dual.M_primal;
        }

        return *this;
    }

    basis_type const& operator()() const
    {
        return M_primal.basis();
    }

    //@}

    /** @name Accessors
     */
    //@{

    primal_space_type const& primalSpace() const
    {
        return M_primal;
    }

    basis_type const& basis() const
    {
        return M_primal.basis();
    }

    //@}

    /** @name  Mutators
     */
    //@{


    //@}

    /** @name  Methods
     */
    //@{



    //@}



protected:

private:

    primal_space_type M_primal;
};
} // Feel

#endif /* __DualBasis_H */
