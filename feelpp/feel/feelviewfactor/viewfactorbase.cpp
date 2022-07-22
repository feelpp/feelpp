/**
 * @file viewfactorbase.cpp
 * @author Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 * @brief 
 * @version 0.1
 * @date 2022-07-22
 * 
 * @copyright Copyright (c) 2022 Feel++ Consortium
 * @copyright copyright (c) 2022 Université de Strasbourg
 * 
 */

#include <feel/feelviewfactor/viewfactorbase.hpp>

namespace Feel {

ViewFactorBase::value_type
ViewFactorBase::devReciprocity( unsigned int i, unsigned int j ) const
{
    return vf_(i,j) - areas_(j) / areas_(i) * vf_(j,i);
}
ViewFactorBase::value_type
ViewFactorBase::maxDevReciprocity() const
{
    return vf_.array().rowwise()
}


} // namespace Feel