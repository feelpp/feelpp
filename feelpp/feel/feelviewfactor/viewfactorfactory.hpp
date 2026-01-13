/**
 * @file viewfactorfactory.hpp
 * @author Christophe Prud'homme (christophe.prudhomme@cemosis.fr)
 * @brief view factory factory
 * @version 0.1
 * @date 2022-07-28
 *
 * @copyright Copyright (c) 2022 Feel++ Consortium
 * @copyright Copyright (c) 2022 Université de Strasbourg
 *
 */
#pragma once

#include <memory>
#include <feel/feelviewfactor/viewfactorbase.hpp>

namespace Feel
{
/**
 * @brief view factor producer factory
 * @ingroup ViewFactor
 */
template<typename MeshType>
class ViewFactorProducerFactory
{
public:
    virtual ~ViewFactorProducerFactory() = default;
    virtual std::unique_ptr<ViewFactorBase<MeshType>> create() = 0;
};

} // namespace Feel
