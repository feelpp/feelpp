// SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
// SPDX-FileCopyrightText: 2026 University of Strasbourg
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <ginac/ginac.h>

namespace Feel::detail
{
/** Return GiNaC's parser table extended with the Feel++ expression language.
 * Stock and bundled providers expose the same function names to expression
 * parsing; the table preserves all provider-supplied functions as well.
 */
GiNaC::prototype_table const& ginacExpressionReader();
} // namespace Feel::detail
