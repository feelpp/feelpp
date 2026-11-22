// SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
// SPDX-FileCopyrightText: 2026 University of Strasbourg
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <functional>
#include <memory>
#include <string>

namespace Feel
{
/** Factory called by Environment after MPI and solver initialization.
 * Environment retains the lease until numerical objects have been cleared,
 * then releases it before solver and MPI finalization.
 */
using EnvironmentRuntimeFactory = std::function<std::shared_ptr<void>()>;

/** Register an optional execution runtime before constructing Environment.
 * @param name Unique runtime name.
 * @param factory Factory returning the Environment-owned runtime lease.
 * Registration and Environment construction occur on the main thread.
 * @throws std::logic_error For duplicate or late registration.
 * @throws std::invalid_argument If name or factory is empty.
 */
void registerEnvironmentRuntime( std::string name, EnvironmentRuntimeFactory factory );
} // namespace Feel
