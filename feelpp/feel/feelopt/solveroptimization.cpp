/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/

#include <feel/feelopt/solveroptimization.hpp>

#if defined( FEELPP_HAS_PETSC_TAO )
#include <feel/feelopt/solveroptimizationpetsc.hpp>
#endif

namespace Feel
{

template<typename T, typename SizeT>
typename SolverOptimization<T, SizeT>::ptrtype
SolverOptimization<T, SizeT>::build(
    std::string const& backend, std::string const& prefix,
    worldcomm_ptr_t const& worldComm, po::variables_map const& vm )
{
    if ( backend != "petsc" )
        throw std::invalid_argument(
            "SolverOptimization::build does not support backend '" + backend + "'" );

#if defined( FEELPP_HAS_PETSC_TAO )
    return std::make_shared<SolverOptimizationPetsc<T, SizeT>>(
        prefix, worldComm, vm );
#else
    static_cast<void>( prefix );
    static_cast<void>( worldComm );
    static_cast<void>( vm );
    throw std::runtime_error(
        "SolverOptimization::build requested PETSc TAO, but this Feel++ build "
        "does not provide it" );
#endif
}

template class SolverOptimization<double>;

} // namespace Feel
