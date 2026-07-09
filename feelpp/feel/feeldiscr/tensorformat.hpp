/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

    SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>

    SPDX-FileCopyrightText: 2026 University of Strasbourg

    SPDX-License-Identifier: LGPL-3.0-or-later
*/
/**
 * @file tensorformat.hpp
 * @brief Symmetric tensor storage, output-order, and scaling utilities.
 *
 * This header centralizes field-level symmetric tensor conventions used by
 * matrix fields and exporters. It deliberately depends only on low-level
 * Feel++ traits for the legacy packed storage index and does not pull in VF
 * expression machinery.
 */
#ifndef FEELPP_FEELDISCR_TENSORFORMAT_HPP
#define FEELPP_FEELDISCR_TENSORFORMAT_HPP 1

#include <array>
#include <concepts>
#include <limits>
#include <stdexcept>

#include <feel/feelpoly/traits.hpp>

namespace Feel
{

/**
 * @brief Output order for packed symmetric tensor components.
 *
 * `Storage` and `DiagonalFirst` are compact orders with 3 components in 2D and
 * 6 components in 3D. The `*Tensor6` values describe six-slot exporter tuples
 * used for both 2D and 3D output, with unused 2D z-components left empty by
 * callers.
 */
enum class SymmetricTensorOrder
{
    /** Legacy packed storage order: 2D `xx,xy,yy`; 3D `xx,xy,xz,yy,yz,zz`. */
    Storage,
    /** Compact diagonal-first order: 2D `xx,yy,xy`; 3D `xx,yy,zz,xy,xz,yz`. */
    DiagonalFirst,
    /** Current VTK tensor-six tuple: `xx,yy,zz,xy,yz,xz`. */
    VtkTensor6,
    /** Current Ensight tensor-six tuple: `xx,yy,zz,xy,xz,yz`. */
    EnsightTensor6,
    /** Current XDMF tensor-six tuple matching 3D storage slots: `xx,xy,xz,yy,yz,zz`. */
    XdmfTensor6
};

/**
 * @brief Scaling convention applied to off-diagonal symmetric tensor entries.
 */
enum class SymmetricTensorScaling
{
    /** No shear scaling; off-diagonal entries keep their tensor value. */
    Tensor,
    /** Engineering shear convention; off-diagonal entries are multiplied by 2. */
    EngineeringShear,
    /** Mandel convention; off-diagonal entries are multiplied by sqrt(2). */
    Mandel
};

/**
 * @brief Field-level symmetric tensor format.
 *
 * The order and scaling are intentionally separate so callers must state both
 * positional order and shear convention. This avoids conflating storage order,
 * exporter tuple order, and Voigt/Mandel-style scaling.
 */
struct SymmetricTensorFormat
{
    /** Requested output order for symmetric tensor components. */
    SymmetricTensorOrder order = SymmetricTensorOrder::Storage;
    /** Requested scaling for off-diagonal components. */
    SymmetricTensorScaling scaling = SymmetricTensorScaling::Tensor;
};

/**
 * @brief Structural concept for objects carrying a symmetric tensor format.
 *
 * Keep this concept local to `tensorformat.hpp` while it only constrains
 * tensor-format helper overloads. Move or mirror it in `feeldiscr/concepts.hpp`
 * only if broader discretization algorithms start depending on it.
 */
template<typename T>
concept SymmetricTensorFormatLike = requires( T const& format )
{
    { format.order } -> std::convertible_to<SymmetricTensorOrder>;
    { format.scaling } -> std::convertible_to<SymmetricTensorScaling>;
};

/**
 * @brief Sentinel value for unused entries in fixed-size mapping arrays.
 *
 * @return The maximum value of `uint16_type`, never a valid output slot.
 */
inline constexpr uint16_type
symmetricTensorInvalidSlot()
{
    return (std::numeric_limits<uint16_type>::max)();
}

/**
 * @brief Check whether a tensor dimension is supported by these helpers.
 *
 * @param n Tensor dimension.
 * @return `true` for dimensions 2 and 3, `false` otherwise.
 */
inline constexpr bool
symmetricTensorSupportedDimension( uint16_type n )
{
    return n == 2 || n == 3;
}

/**
 * @brief Check whether an order describes a six-component exporter tuple.
 *
 * @param order Symmetric tensor order.
 * @return `true` for VTK, Ensight, and XDMF tensor-six orders.
 */
inline constexpr bool
symmetricTensorIsTensor6Order( SymmetricTensorOrder order )
{
    return order == SymmetricTensorOrder::VtkTensor6 ||
           order == SymmetricTensorOrder::EnsightTensor6 ||
           order == SymmetricTensorOrder::XdmfTensor6;
}

/**
 * @brief Stable metadata name for a symmetric tensor output order.
 *
 * @param order Symmetric tensor order.
 * @return String literal matching the enum value name.
 * @throws std::invalid_argument if `order` is unsupported.
 */
inline constexpr char const*
symmetricTensorOrderName( SymmetricTensorOrder order )
{
    return order == SymmetricTensorOrder::Storage ? "Storage" :
           order == SymmetricTensorOrder::DiagonalFirst ? "DiagonalFirst" :
           order == SymmetricTensorOrder::VtkTensor6 ? "VtkTensor6" :
           order == SymmetricTensorOrder::EnsightTensor6 ? "EnsightTensor6" :
           order == SymmetricTensorOrder::XdmfTensor6 ? "XdmfTensor6" :
           throw std::invalid_argument( "unsupported symmetric tensor order" );
}

/**
 * @brief Stable metadata name for a symmetric tensor scaling convention.
 *
 * @param scaling Symmetric tensor scaling convention.
 * @return String literal matching the enum value name.
 * @throws std::invalid_argument if `scaling` is unsupported.
 */
inline constexpr char const*
symmetricTensorScalingName( SymmetricTensorScaling scaling )
{
    return scaling == SymmetricTensorScaling::Tensor ? "Tensor" :
           scaling == SymmetricTensorScaling::EngineeringShear ? "EngineeringShear" :
           scaling == SymmetricTensorScaling::Mandel ? "Mandel" :
           throw std::invalid_argument( "unsupported symmetric tensor scaling" );
}

/**
 * @brief Number of compact packed-storage components for a symmetric tensor.
 *
 * @param n Tensor dimension, currently 2 or 3.
 * @return 3 in 2D and 6 in 3D.
 * @throws std::invalid_argument if `n` is unsupported.
 */
inline constexpr uint16_type
symmetricTensorStorageComponentCount( uint16_type n )
{
    return n == 2 ? 3 :
           n == 3 ? 6 :
           throw std::invalid_argument( "unsupported symmetric tensor dimension" );
}

/**
 * @brief Number of output components for a dimension/order pair.
 *
 * Compact orders return 3 components in 2D and 6 in 3D. Tensor-six exporter
 * orders always return 6 components.
 *
 * @param n Tensor dimension, currently 2 or 3.
 * @param order Requested component order.
 * @return Number of output components.
 * @throws std::invalid_argument if `n` is unsupported.
 */
inline constexpr uint16_type
symmetricTensorComponentCount( uint16_type n, SymmetricTensorOrder order )
{
    return !symmetricTensorSupportedDimension( n ) ?
               throw std::invalid_argument( "unsupported symmetric tensor dimension" ) :
           symmetricTensorIsTensor6Order( order ) ? 6 :
           symmetricTensorStorageComponentCount( n );
}

/**
 * @brief Number of output components for a dimension/format pair.
 *
 * @tparam FormatT Type satisfying `SymmetricTensorFormatLike`.
 * @param n Tensor dimension, currently 2 or 3.
 * @param format Tensor format carrying order and scaling.
 * @return Number of output components implied by `format.order`.
 * @throws std::invalid_argument if `n` is unsupported.
 */
template<SymmetricTensorFormatLike FormatT>
inline constexpr uint16_type
symmetricTensorComponentCount( uint16_type n, FormatT const& format )
{
    return symmetricTensorComponentCount( n, format.order );
}

/**
 * @brief Legacy compact packed-storage slot for a symmetric tensor component.
 *
 * This is the public field-level wrapper around `Feel::detail::symmetricIndex`.
 * It preserves the current storage order: 2D `xx,xy,yy` and 3D
 * `xx,xy,xz,yy,yz,zz`.
 *
 * @param i Row component index.
 * @param j Column component index.
 * @param n Tensor dimension, currently 2 or 3.
 * @return Compact packed-storage slot.
 * @throws std::invalid_argument if `n`, `i`, or `j` is out of range.
 */
inline constexpr uint16_type
symmetricTensorStorageIndex( uint16_type i, uint16_type j, uint16_type n )
{
    return !symmetricTensorSupportedDimension( n ) ?
               throw std::invalid_argument( "unsupported symmetric tensor dimension" ) :
           ( i >= n || j >= n ) ?
               throw std::invalid_argument( "symmetric tensor component is out of range" ) :
           Feel::detail::symmetricIndex( i, j, n );
}

/**
 * @brief Component indices represented by a compact packed-storage slot.
 *
 * @param storageSlot Compact storage slot.
 * @param n Tensor dimension, currently 2 or 3.
 * @return `{i,j}` with `i <= j`.
 * @throws std::invalid_argument if `n` or `storageSlot` is out of range.
 */
inline constexpr std::array<uint16_type,2>
symmetricTensorStorageComponent( uint16_type storageSlot, uint16_type n )
{
    return !symmetricTensorSupportedDimension( n ) ?
               throw std::invalid_argument( "unsupported symmetric tensor dimension" ) :
           storageSlot >= symmetricTensorStorageComponentCount( n ) ?
               throw std::invalid_argument( "symmetric tensor storage slot is out of range" ) :
           n == 2 ?
               ( storageSlot == 0 ? std::array<uint16_type,2>{{0,0}} :
                 storageSlot == 1 ? std::array<uint16_type,2>{{0,1}} :
                                    std::array<uint16_type,2>{{1,1}} ) :
               ( storageSlot == 0 ? std::array<uint16_type,2>{{0,0}} :
                 storageSlot == 1 ? std::array<uint16_type,2>{{0,1}} :
                 storageSlot == 2 ? std::array<uint16_type,2>{{0,2}} :
                 storageSlot == 3 ? std::array<uint16_type,2>{{1,1}} :
                 storageSlot == 4 ? std::array<uint16_type,2>{{1,2}} :
                                    std::array<uint16_type,2>{{2,2}} );
}

/**
 * @brief Check whether a tensor component is off-diagonal.
 *
 * @param i Row component index.
 * @param j Column component index.
 * @param n Tensor dimension, currently 2 or 3.
 * @return `true` when `i != j`.
 * @throws std::invalid_argument if `n`, `i`, or `j` is out of range.
 */
inline constexpr bool
symmetricTensorIsShearComponent( uint16_type i, uint16_type j, uint16_type n )
{
    return !symmetricTensorSupportedDimension( n ) ?
               throw std::invalid_argument( "unsupported symmetric tensor dimension" ) :
           ( i >= n || j >= n ) ?
               throw std::invalid_argument( "symmetric tensor component is out of range" ) :
           i != j;
}

/**
 * @brief Check whether a compact storage slot represents an off-diagonal component.
 *
 * @param storageSlot Compact storage slot.
 * @param n Tensor dimension, currently 2 or 3.
 * @return `true` when the slot maps to an off-diagonal component.
 * @throws std::invalid_argument if `n` or `storageSlot` is out of range.
 */
inline constexpr bool
symmetricTensorIsShearStorageSlot( uint16_type storageSlot, uint16_type n )
{
    auto const ij = symmetricTensorStorageComponent( storageSlot, n );
    return ij[0] != ij[1];
}

/**
 * @brief Map a compact storage slot to an output slot.
 *
 * For compact output orders, the output slot is in a compact vector. For
 * tensor-six exporter orders, the output slot is in the six-component tuple.
 *
 * @param storageSlot Compact storage slot in the legacy field storage order.
 * @param n Tensor dimension, currently 2 or 3.
 * @param order Requested output order.
 * @return Output slot for `storageSlot`.
 * @throws std::invalid_argument if `n` or `storageSlot` is out of range.
 */
inline constexpr uint16_type
symmetricTensorOutputSlotFromStorage( uint16_type storageSlot, uint16_type n, SymmetricTensorOrder order )
{
    return storageSlot >= symmetricTensorStorageComponentCount( n ) ?
               throw std::invalid_argument( "symmetric tensor storage slot is out of range" ) :
           order == SymmetricTensorOrder::Storage ?
               storageSlot :
           order == SymmetricTensorOrder::DiagonalFirst ?
               ( n == 2 ?
                 ( storageSlot == 0 ? 0 :
                   storageSlot == 1 ? 2 : 1 ) :
                 ( storageSlot == 0 ? 0 :
                   storageSlot == 1 ? 3 :
                   storageSlot == 2 ? 4 :
                   storageSlot == 3 ? 1 :
                   storageSlot == 4 ? 5 : 2 ) ) :
           order == SymmetricTensorOrder::VtkTensor6 ?
               ( n == 2 ?
                 ( storageSlot == 0 ? 0 :
                   storageSlot == 1 ? 3 : 1 ) :
                 ( storageSlot == 0 ? 0 :
                   storageSlot == 1 ? 3 :
                   storageSlot == 2 ? 5 :
                   storageSlot == 3 ? 1 :
                   storageSlot == 4 ? 4 : 2 ) ) :
           order == SymmetricTensorOrder::EnsightTensor6 ?
               ( n == 2 ?
                 ( storageSlot == 0 ? 0 :
                   storageSlot == 1 ? 3 : 1 ) :
                 ( storageSlot == 0 ? 0 :
                   storageSlot == 1 ? 3 :
                   storageSlot == 2 ? 4 :
                   storageSlot == 3 ? 1 :
                   storageSlot == 4 ? 5 : 2 ) ) :
           order == SymmetricTensorOrder::XdmfTensor6 ?
               ( n == 2 ?
                 ( storageSlot == 0 ? 0 :
                   storageSlot == 1 ? 1 : 3 ) :
                 storageSlot ) :
           throw std::invalid_argument( "unsupported symmetric tensor order" );
}

/**
 * @brief Map a compact storage slot to an output slot for a format object.
 *
 * @tparam FormatT Type satisfying `SymmetricTensorFormatLike`.
 * @param storageSlot Compact storage slot in the legacy field storage order.
 * @param n Tensor dimension, currently 2 or 3.
 * @param format Tensor format carrying order and scaling.
 * @return Output slot for `storageSlot` using `format.order`.
 * @throws std::invalid_argument if `n` or `storageSlot` is out of range.
 */
template<SymmetricTensorFormatLike FormatT>
inline constexpr uint16_type
symmetricTensorOutputSlotFromStorage( uint16_type storageSlot, uint16_type n, FormatT const& format )
{
    return symmetricTensorOutputSlotFromStorage( storageSlot, n, format.order );
}

/**
 * @brief Map a symmetric tensor component to an output slot.
 *
 * @param i Row component index.
 * @param j Column component index.
 * @param n Tensor dimension, currently 2 or 3.
 * @param order Requested output order.
 * @return Output slot for component `(i,j)`.
 * @throws std::invalid_argument if `n`, `i`, or `j` is out of range.
 */
inline constexpr uint16_type
symmetricTensorOutputSlot( uint16_type i, uint16_type j, uint16_type n, SymmetricTensorOrder order )
{
    return symmetricTensorOutputSlotFromStorage( symmetricTensorStorageIndex( i, j, n ), n, order );
}

/**
 * @brief Map a symmetric tensor component to an output slot for a format object.
 *
 * @tparam FormatT Type satisfying `SymmetricTensorFormatLike`.
 * @param i Row component index.
 * @param j Column component index.
 * @param n Tensor dimension, currently 2 or 3.
 * @param format Tensor format carrying order and scaling.
 * @return Output slot for component `(i,j)` using `format.order`.
 * @throws std::invalid_argument if `n`, `i`, or `j` is out of range.
 */
template<SymmetricTensorFormatLike FormatT>
inline constexpr uint16_type
symmetricTensorOutputSlot( uint16_type i, uint16_type j, uint16_type n, FormatT const& format )
{
    return symmetricTensorOutputSlot( i, j, n, format.order );
}

/**
 * @brief Fixed-size compact-storage-to-output-slot map.
 *
 * The first `symmetricTensorStorageComponentCount(n)` entries are valid. Any
 * remaining entries are filled with `symmetricTensorInvalidSlot()`.
 *
 * @param n Tensor dimension, currently 2 or 3.
 * @param order Requested output order.
 * @return Six-entry mapping array from compact storage slot to output slot.
 * @throws std::invalid_argument if `n` is unsupported.
 */
inline constexpr std::array<uint16_type,6>
symmetricTensorStorageToOutputMap( uint16_type n, SymmetricTensorOrder order )
{
    return n == 2 ?
           std::array<uint16_type,6>{{
               symmetricTensorOutputSlotFromStorage( 0, n, order ),
               symmetricTensorOutputSlotFromStorage( 1, n, order ),
               symmetricTensorOutputSlotFromStorage( 2, n, order ),
               symmetricTensorInvalidSlot(),
               symmetricTensorInvalidSlot(),
               symmetricTensorInvalidSlot()
           }} :
           n == 3 ?
           std::array<uint16_type,6>{{
               symmetricTensorOutputSlotFromStorage( 0, n, order ),
               symmetricTensorOutputSlotFromStorage( 1, n, order ),
               symmetricTensorOutputSlotFromStorage( 2, n, order ),
               symmetricTensorOutputSlotFromStorage( 3, n, order ),
               symmetricTensorOutputSlotFromStorage( 4, n, order ),
               symmetricTensorOutputSlotFromStorage( 5, n, order )
           }} :
           throw std::invalid_argument( "unsupported symmetric tensor dimension" );
}

/**
 * @brief Component label for an output slot.
 *
 * @param outputSlot Slot in the requested output order.
 * @param n Tensor dimension, currently 2 or 3.
 * @param order Requested output order.
 * @return Static label such as `"xx"`, `"xy"`, or `"zz"`.
 * @throws std::invalid_argument if `n` or `outputSlot` is out of range.
 */
inline constexpr char const*
symmetricTensorComponentLabel( uint16_type outputSlot, uint16_type n, SymmetricTensorOrder order )
{
    return outputSlot >= symmetricTensorComponentCount( n, order ) ?
               throw std::invalid_argument( "symmetric tensor output slot is out of range" ) :
           order == SymmetricTensorOrder::Storage ?
               ( n == 2 ?
                 ( outputSlot == 0 ? "xx" :
                   outputSlot == 1 ? "xy" : "yy" ) :
                 ( outputSlot == 0 ? "xx" :
                   outputSlot == 1 ? "xy" :
                   outputSlot == 2 ? "xz" :
                   outputSlot == 3 ? "yy" :
                   outputSlot == 4 ? "yz" : "zz" ) ) :
           order == SymmetricTensorOrder::DiagonalFirst ?
               ( n == 2 ?
                 ( outputSlot == 0 ? "xx" :
                   outputSlot == 1 ? "yy" : "xy" ) :
                 ( outputSlot == 0 ? "xx" :
                   outputSlot == 1 ? "yy" :
                   outputSlot == 2 ? "zz" :
                   outputSlot == 3 ? "xy" :
                   outputSlot == 4 ? "xz" : "yz" ) ) :
           order == SymmetricTensorOrder::VtkTensor6 ?
               ( outputSlot == 0 ? "xx" :
                 outputSlot == 1 ? "yy" :
                 outputSlot == 2 ? "zz" :
                 outputSlot == 3 ? "xy" :
                 outputSlot == 4 ? "yz" : "xz" ) :
           order == SymmetricTensorOrder::EnsightTensor6 ?
               ( outputSlot == 0 ? "xx" :
                 outputSlot == 1 ? "yy" :
                 outputSlot == 2 ? "zz" :
                 outputSlot == 3 ? "xy" :
                 outputSlot == 4 ? "xz" : "yz" ) :
           order == SymmetricTensorOrder::XdmfTensor6 ?
               ( outputSlot == 0 ? "xx" :
                 outputSlot == 1 ? "xy" :
                 outputSlot == 2 ? "xz" :
                 outputSlot == 3 ? "yy" :
                 outputSlot == 4 ? "yz" : "zz" ) :
           throw std::invalid_argument( "unsupported symmetric tensor order" );
}

/**
 * @brief Component label for an output slot in a format object.
 *
 * @tparam FormatT Type satisfying `SymmetricTensorFormatLike`.
 * @param outputSlot Slot in the requested output order.
 * @param n Tensor dimension, currently 2 or 3.
 * @param format Tensor format carrying order and scaling.
 * @return Static label such as `"xx"`, `"xy"`, or `"zz"`.
 * @throws std::invalid_argument if `n` or `outputSlot` is out of range.
 */
template<SymmetricTensorFormatLike FormatT>
inline constexpr char const*
symmetricTensorComponentLabel( uint16_type outputSlot, uint16_type n, FormatT const& format )
{
    return symmetricTensorComponentLabel( outputSlot, n, format.order );
}

/**
 * @brief Scaling factor for a diagonal or off-diagonal component.
 *
 * @tparam T Floating-point-like return type.
 * @param isShear Whether the component is off-diagonal.
 * @param scaling Requested shear scaling convention.
 * @return `1` for tensor scaling and all diagonal components, `2` for
 * engineering shear, and `sqrt(2)` for Mandel shear.
 * @throws std::invalid_argument if `scaling` is unsupported.
 */
template<typename T = double>
inline constexpr T
symmetricTensorScale( bool isShear, SymmetricTensorScaling scaling )
{
    return !isShear ? T( 1 ) :
           scaling == SymmetricTensorScaling::Tensor ? T( 1 ) :
           scaling == SymmetricTensorScaling::EngineeringShear ? T( 2 ) :
           scaling == SymmetricTensorScaling::Mandel ? static_cast<T>( 1.41421356237309504880168872420969808L ) :
           throw std::invalid_argument( "unsupported symmetric tensor scaling" );
}

/**
 * @brief Scaling factor for a tensor component.
 *
 * @tparam T Floating-point-like return type.
 * @param i Row component index.
 * @param j Column component index.
 * @param n Tensor dimension, currently 2 or 3.
 * @param scaling Requested shear scaling convention.
 * @return Scaling factor for component `(i,j)`.
 * @throws std::invalid_argument if `n`, `i`, or `j` is out of range.
 */
template<typename T = double>
inline constexpr T
symmetricTensorComponentScale( uint16_type i, uint16_type j, uint16_type n, SymmetricTensorScaling scaling )
{
    return symmetricTensorScale<T>( symmetricTensorIsShearComponent( i, j, n ), scaling );
}

/**
 * @brief Scaling factor for a compact storage slot.
 *
 * @tparam T Floating-point-like return type.
 * @param storageSlot Compact storage slot in the legacy field storage order.
 * @param n Tensor dimension, currently 2 or 3.
 * @param scaling Requested shear scaling convention.
 * @return Scaling factor for the component represented by `storageSlot`.
 * @throws std::invalid_argument if `n` or `storageSlot` is out of range.
 */
template<typename T = double>
inline constexpr T
symmetricTensorStorageScale( uint16_type storageSlot, uint16_type n, SymmetricTensorScaling scaling )
{
    return symmetricTensorScale<T>( symmetricTensorIsShearStorageSlot( storageSlot, n ), scaling );
}

} // namespace Feel

#endif // FEELPP_FEELDISCR_TENSORFORMAT_HPP
