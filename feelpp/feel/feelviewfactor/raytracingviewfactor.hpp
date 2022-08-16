/**
 * @file raytracingviewfactor.hpp
 * @author Christophe Prud'homme (christophe.prudhomme@cemosis.fr)
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

template <typename MeshType>
class RayTracingViewFactor : public ViewFactorBase<MeshType>
{
public:
    using value_type = double;
    RayTracingViewFactor() = default;
    RayTracingViewFactor( std::vector<std::string> const& list_of_bdys ) 
        : 
        ViewFactorBase( list_of_bdys )
        {}
    RayTracingViewFactor( const RayTracingViewFactor& ) = default;
    RayTracingViewFactor( RayTracingViewFactor&& ) = default;
    RayTracingViewFactor& operator=( const RayTracingViewFactor& ) = default;
    RayTracingViewFactor& operator=( RayTracingViewFactor&& ) = default;
    ~RayTracingViewFactor() = default;
    void init( std::vector<std::string> const& list_of_bdys ) { ViewFactorBase::init( list_of_bdys ); }
};

} // namespace Feel
};
} // namespace Feel