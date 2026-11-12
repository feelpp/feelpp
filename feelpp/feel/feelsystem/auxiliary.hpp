/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@cemosis.fr>
 *
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 *
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */

#ifndef FEELPP_FEELSYSTEM_AUXILIARY_HPP
#define FEELPP_FEELSYSTEM_AUXILIARY_HPP 1

/**
 * @file auxiliary.hpp
 * @brief Backend-independent auxiliary operators assembled by FeelSystem.
 */

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <feel/feelalg/backend.hpp>
#include <feel/feelalg/preconditioner.hpp>

namespace Feel::systems
{

/** @brief Mathematical meaning of an auxiliary operator. */
enum class AuxiliaryOperatorKind
{
    local_neumann
};

/** @brief Terms used to assemble an auxiliary operator. */
enum class AuxiliaryOperatorSource
{
    explicit_terms,
    lhs_terms,
    inherited_problem
};

/**
 * @brief Backend-neutral declaration of a named auxiliary operator.
 *
 * The descriptor deliberately carries no PETSc or preconditioner options. It
 * says what to assemble; a later backend adapter decides how to consume it.
 */
struct AuxiliaryOperatorDescriptor
{
    std::string name;
    AuxiliaryOperatorKind kind = AuxiliaryOperatorKind::local_neumann;
    AuxiliaryOperatorSource source = AuxiliaryOperatorSource::explicit_terms;
};

/**
 * @brief Per-rank unassembled auxiliary matrix and its global overlap mapping.
 *
 * `matrixPtr()` is sequential on every MPI rank. Entry `i` of
 * `globalDofIds()` identifies local row/column `i` in the distributed system's
 * global numbering. This is the representation required by overlapping-domain
 * preconditioners while remaining independent of a particular backend API.
 */
template <typename ValueType>
class AssembledAuxiliaryOperator
{
  public:
    using value_type = ValueType;
    using backend_type = Backend<value_type>;
    using backend_ptrtype = std::shared_ptr<backend_type>;
    using matrix_ptrtype = typename backend_type::sparse_matrix_ptrtype;
    using size_type = Feel::size_type;

    AssembledAuxiliaryOperator() = default;

    AssembledAuxiliaryOperator( AuxiliaryOperatorDescriptor descriptor,
                                backend_ptrtype backend,
                                matrix_ptrtype matrix,
                                std::vector<size_type> globalDofIds )
        : M_descriptor( std::move( descriptor ) ),
          M_backend( std::move( backend ) ),
          M_matrix( std::move( matrix ) ),
          M_globalDofIds( std::move( globalDofIds ) )
    {
    }

    AuxiliaryOperatorDescriptor const& descriptor() const noexcept { return M_descriptor; }
    std::string const& name() const noexcept { return M_descriptor.name; }
    AuxiliaryOperatorKind kind() const noexcept { return M_descriptor.kind; }
    backend_ptrtype const& backend() const noexcept { return M_backend; }
    matrix_ptrtype const& matrixPtr() const noexcept { return M_matrix; }
    std::vector<size_type> const& globalDofIds() const noexcept { return M_globalDofIds; }
    std::size_t localSize() const noexcept { return M_globalDofIds.size(); }

  private:
    AuxiliaryOperatorDescriptor M_descriptor;
    backend_ptrtype M_backend;
    matrix_ptrtype M_matrix;
    std::vector<size_type> M_globalDofIds;
};

namespace detail
{

struct AuxiliaryDataMap
{
    std::shared_ptr<DataMap<>> localMap;
    std::vector<size_type> globalDofIds;
};

/**
 * @brief Build a sequential tagged map matching a distributed process-local map.
 */
inline AuxiliaryDataMap
makeAuxiliaryDataMap( std::shared_ptr<DataMap<>> const& distributedMap )
{
    CHECK( distributedMap ) << "an auxiliary operator requires a source data map";
    auto const localSize = distributedMap->nLocalDofWithGhost();
    auto localMap = std::make_shared<DataMap<>>( localSize, localSize, Environment::worldCommSeqPtr() );
    std::vector<size_type> globalDofIds( localSize );

    int const nTags = distributedMap->numberOfDofIdToContainerId();
    if ( nTags == 0 )
    {
        for ( size_type i = 0; i < localSize; ++i )
            globalDofIds[i] = distributedMap->mapGlobalProcessToGlobalCluster( i );
        return { std::move( localMap ), std::move( globalDofIds ) };
    }

    localMap->initNumberOfDofIdToContainerId( nTags );
    size_type localIndex = 0;
    for ( int tag = 0; tag < nTags; ++tag )
    {
        auto const& distributedIds = distributedMap->dofIdToContainerId( tag );
        localMap->initDofIdToContainerId( tag, distributedIds.size() );
        auto& localIds = localMap->dofIdToContainerIdRef( tag );
        for ( size_type dof = 0; dof < distributedIds.size(); ++dof )
        {
            CHECK_LT( localIndex, localSize ) << "inconsistent tagged auxiliary data map";
            localIds[dof] = localIndex;
            globalDofIds[localIndex] = distributedMap->mapGlobalProcessToGlobalCluster( distributedIds[dof] );
            ++localIndex;
        }
    }
    CHECK_EQ( localIndex, localSize ) << "tagged source map does not cover all local degrees of freedom";
    return { std::move( localMap ), std::move( globalDofIds ) };
}

/**
 * @brief Attach semantic FeelSystem auxiliary operators to a preconditioner.
 *
 * HPDDM consumes one floating local Neumann operator. Other preconditioners
 * ignore these operators, so merely declaring one does not change the solver
 * unless HPDDM is selected.
 */
template <typename ValueType>
void
attachAuxiliaryOperators(
    std::shared_ptr<Preconditioner<ValueType>> const& preconditioner,
    std::vector<AssembledAuxiliaryOperator<ValueType>> const& operators )
{
    CHECK( preconditioner ) << "auxiliary operators require a preconditioner";
    if ( preconditioner->type() != HPDDM_PRECOND )
        return;

    AssembledAuxiliaryOperator<ValueType> const* localNeumann = nullptr;
    for ( auto const& op : operators )
    {
        if ( op.kind() != AuxiliaryOperatorKind::local_neumann )
            continue;
        CHECK( localNeumann == nullptr )
            << "HPDDM accepts exactly one local Neumann operator per preconditioner";
        localNeumann = &op;
    }
    CHECK( localNeumann != nullptr )
        << "HPDDM auxiliary attachment requires a local Neumann operator";
    preconditioner->attachLocalNeumannMatrix( localNeumann->matrixPtr(),
                                              localNeumann->globalDofIds() );
}

/** @brief Build the explicit preconditioner needed for auxiliary attachment. */
template <typename ValueType>
std::shared_ptr<Preconditioner<ValueType>>
makeAuxiliaryPreconditioner(
    std::shared_ptr<Backend<ValueType>> const& backend,
    typename Backend<ValueType>::sparse_matrix_ptrtype const& matrix,
    std::string const& pc,
    std::string const& factorPackage )
{
    CHECK( backend ) << "auxiliary operators require a backend";
    auto preconditioner = Feel::preconditioner(
        _prefix = backend->prefix(),
        _matrix = matrix,
        _pc = pcTypeConvertStrToEnum( pc ),
        _pcfactormatsolverpackage = matSolverPackageConvertStrToEnum( factorPackage ),
        _backend = backend,
        _worldcomm = backend->worldCommPtr() );
    return preconditioner;
}

} // namespace detail

} // namespace Feel::systems

#endif // FEELPP_FEELSYSTEM_AUXILIARY_HPP
