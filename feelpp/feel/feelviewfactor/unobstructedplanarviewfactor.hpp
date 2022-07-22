/**
 * @file unobstructedplanarviewfactor.hpp
 * @author Christophe Prud'homme
 * @brief
 * @version 0.1
 * @date 2022-07-22
 *
 * @copyright Copyright (c) 2022 Feel++ Consortium
 * @copyright Copyright (c) 2022 Université de Strasbourg
 *
 */
#pragma once

#include <feel/feelviewfactor/viewfactorbase.hpp>

namespace Feel {

/**
 * @brief  Computes the view factors for planar faces in unobstructed radiative heat transfer
 *
 */
class UnobstructedPlanarViewFactor : public ViewFactorBase
{
public:
    using value_type = double;
    UnobstructedPlanarViewFactor() = default;
    UnobstructedPlanarViewFactor( std::vector<std::string> const& list_of_bdys ) 
        : 
        ViewFactorBase( list_of_bdys )
        {}
    UnobstructedPlanarViewFactor( const UnobstructedPlanarViewFactor& ) = default;
    UnobstructedPlanarViewFactor( UnobstructedPlanarViewFactor&& ) = default;
    UnobstructedPlanarViewFactor& operator=( const UnobstructedPlanarViewFactor& ) = default;
    UnobstructedPlanarViewFactor& operator=( UnobstructedPlanarViewFactor&& ) = default;
    ~UnobstructedPlanarViewFactor() = default;
    void init( std::vector<std::string> const& list_of_bdys ) { ViewFactorBase::init( list_of_bdys ); }
};

} // namespace Feel