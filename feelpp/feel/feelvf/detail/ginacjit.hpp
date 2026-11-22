// SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
// SPDX-FileCopyrightText: 2026 University of Strasbourg
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <ginac/ginac.h>
#include <string>

namespace boost::mpi
{
class communicator;
}

/** Feel++ expression compilation and shared-filesystem MPI coordination.
 * GiNaC supplies symbolic expressions only. Compiled modules and their cache
 * lifetime belong to this adapter, independently of the GiNaC provider.
 */
namespace Feel::detail::ginacjit
{
using GiNaC::ex;
using GiNaC::lst;
using GiNaC::symbol;

using GiNaC::FUNCP_1P;
using GiNaC::FUNCP_2P;
using GiNaC::FUNCP_CUBA;

/** Compile a scalar expression collectively on the default communicator. */
void compile_ex( ex const& expression, symbol const& parameter, FUNCP_1P& function, std::string filename = "" );

/** Compile on communicator rank zero and link on every rank; files must be shared. */
void compile_ex( ex const& expression, symbol const& parameter, FUNCP_1P& function, std::string filename, boost::mpi::communicator const& communicator );

/** Compile a two-argument scalar expression on the default communicator. */
void compile_ex( ex const& expression, symbol const& firstParameter, symbol const& secondParameter, FUNCP_2P& function, std::string filename = "" );

/** Compile a two-argument scalar expression collectively on the supplied communicator. */
void compile_ex( ex const& expression, symbol const& firstParameter, symbol const& secondParameter, FUNCP_2P& function, std::string filename, boost::mpi::communicator const& communicator );

/** Compile multiple outputs, optionally using an already content-addressed filename. */
void compile_ex( lst const& expressions, lst const& parameters, FUNCP_CUBA& function, std::string filename = "", bool contentAddressed = false );

/** Compile multiple outputs collectively; all ranks supply identical expressions and paths. */
void compile_ex( lst const& expressions, lst const& parameters, FUNCP_CUBA& function, std::string filename, bool contentAddressed, boost::mpi::communicator const& communicator );

/** Link a previously compiled one-argument module. */
void link_ex( std::string filename, FUNCP_1P& function );

/** Link a previously compiled two-argument module. */
void link_ex( std::string filename, FUNCP_2P& function );

/** Link a previously compiled multi-output module. */
void link_ex( std::string filename, FUNCP_CUBA& function );

/** Release modules loaded from this filename. */
void unlink_ex( std::string filename );

/** Release adapter-owned modules before runtime teardown. */
void cleanup_ex( bool verbose );
} // namespace Feel::detail::ginacjit
