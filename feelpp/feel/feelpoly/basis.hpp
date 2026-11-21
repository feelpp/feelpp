/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2006 EPFL
 * SPDX-FileCopyrightText: 2007 Université Joseph Fourier Grenoble 1
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
/** @file basis.hpp
 * @brief Lazy immutable differentiation operators for reachable legacy descriptors.
 */
#ifndef FEELPP_FEELPOLY_BASIS_HPP
#define FEELPP_FEELPOLY_BASIS_HPP 1

#include <atomic>
#include <feel/feelpoly/referenceoperators.hpp>
#include <span>

namespace Feel
{
/** @brief Preserve static derivative access without evaluating a partly constructed primal.
 * @tparam Tag Complete mathematical configuration identity of the legacy primal.
 * @tparam T Scalar type.
 * @note Construction registers a stateless typed provider only. First derivative access
 * constructs a fully initialized default primal, then publishes all directions through
 * thread-safe function-local initialization. Failed construction is retried.
 */
template <typename Tag, typename T>
class Basis
{
  public:
    /** Scalar type of the legacy facade. */
    using value_type = T;
    /** Existing compatibility return representation. */
    using matrix_type = ublas::matrix<T, ublas::row_major>;

    /** @brief Register a provider without reading uninitialized derived members. */
    template <typename PrimalBasis>
    explicit Basis( PrimalBasis const& primal ) { initDerivation( primal ); }
    /** @brief Copies share the same immutable mathematical configuration. */
    Basis( Basis const& ) = default;
    /** @brief Assignment preserves mathematical configuration identity. */
    Basis& operator=( Basis const& ) = default;
    /** @brief Release the descriptor; operator ownership has independent lifetime. */
    virtual ~Basis() = default;

    /** @brief Borrow a checked direction after constructing a matching primal descriptor. */
    [[nodiscard]] static matrix_type const& d( uint16_type direction )
    {
        auto provider = providerSlot().load( std::memory_order_acquire );
        if ( !provider )
            throw std::logic_error( "construct a matching primal before requesting legacy operators" );
        auto const matrices = provider();
        if ( direction >= matrices.size() )
            throw std::out_of_range( "legacy derivative direction is outside the reference dimension" );
        return matrices[direction];
    }
    /** @brief Established spelling for the same immutable derivative operator. */
    [[nodiscard]] static matrix_type const& derivate( uint16_type direction ) { return d( direction ); }

  protected:
    /** @brief Bind the tag to its primal without eagerly constructing operators. */
    template <typename PrimalBasis>
    static void initDerivation( PrimalBasis const& )
    {
        Provider expected = nullptr;
        providerSlot().compare_exchange_strong( expected, &operatorsFor<PrimalBasis>,
                                                std::memory_order_release, std::memory_order_relaxed );
    }

  private:
    /** Type-erased view retains fixed std::array ownership in the typed provider. */
    using Provider = std::span<matrix_type const> ( * )();

    /** @brief Retain a provider pointer; its local static owns the immutable matrices. */
    static std::atomic<Provider>& providerSlot()
    {
        static std::atomic<Provider> provider{ nullptr };
        return provider;
    }

    /** @brief Construct a complete primal before evaluating its interpolation tables. */
    template <typename PrimalBasis>
    static std::span<matrix_type const> operatorsFor()
    {
        static auto const operators = []
        {
            PrimalBasis primal;
            using Traits = typename PrimalBasis::traits_type;
            typename Traits::diff_pointset_type points( 1 );
            auto values = primal.evaluate( points.points() );
            auto gradients = primal.derivate( points.points() );
            return ::Feel::detail::interpolateReferenceOperators<Traits::nDim, matrix_type>(
                RuntimeOrder{ Traits::nOrder }, values, gradients,
                type_traits<T>::epsilon(), "legacy basis differentiation" );
        }();
        return operators.matrices;
    }
};
} // namespace Feel
#endif // FEELPP_FEELPOLY_BASIS_HPP
