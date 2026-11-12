/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 *
 * SPDX-FileCopyrightText: 2012-2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

/**
   \file systemexplicit.hpp
   \brief Legacy explicit system compatibility wrapper.
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-01-04
 */
#ifndef __SystemExplicit_H
#define __SystemExplicit_H 1

#include <feel/feelsystem/systemlegacy.hpp>

namespace Feel
{
/**
 * \class SystemExplicit
 * \brief describes an explicit system
 *
 * @author Christophe Prud'homme
 * @see Feel::System
 */
template<typename SpaceType>
class SystemExplicit : public System<SpaceType>
{
    /** @brief Base legacy system type. */
    typedef System<SpaceType> super;
public:


    /** @name Constants
     */
    //@{

    /** @brief Topological dimension inherited from the function space. */
    static inline const uint16_type Dim = super::Dim;

    //@}

    /** @name Typedefs
     */
    //@{

    /** @brief Self type. */
    typedef SystemExplicit<SpaceType> system_type;

    /** @brief Scalar value type. */
    typedef typename super::value_type value_type;
    /** @brief Function-space type. */
    typedef typename super::functionspace_type functionspace_type;
    /** @brief Legacy function-space argument type. */
    typedef typename super::functionspace_type functionspace_ptrtype;
    /** @brief Finite-element field type. */
    typedef typename super::element_type element_type;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    /** @brief Construct from a function space and variables map. */
    SystemExplicit( functionspace_ptrtype Xh, po::variables_map const& vm ) : super( Xh, vm ) {}
    /** @brief Copy constructor. */
    SystemExplicit( SystemExplicit const & se ) : super( se ) {}
    /** @brief Destructor. */
    ~SystemExplicit() {}

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}

    /** @name Accessors
     */
    //@{


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

};
} // Feel
#endif /* __SystemExplicit_H */
