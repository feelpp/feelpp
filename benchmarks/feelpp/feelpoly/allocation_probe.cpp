/*
 * SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
 * SPDX-FileCopyrightText: 2026 University of Strasbourg
 * SPDX-License-Identifier: LGPL-3.0-or-later
 */
#include "allocation_probe.hpp"

#include <atomic>
#include <cstddef>
#include <cstdlib>
#include <new>

namespace Feel::benchmarking::allocation_probe
{
/** Whether allocation recording is currently enabled. */
std::atomic<bool> enabled{ false };
/** Cumulative allocation-call counter. */
std::atomic<std::uint64_t> calls{ 0 };
/** Cumulative requested-byte counter. */
std::atomic<std::uint64_t> bytes{ 0 };

/**
 * @brief Record one C++ allocation request when probing is enabled.
 * @param size requested allocation size in bytes
 */
inline void record( std::size_t size ) noexcept
{
    if ( !enabled.load( std::memory_order_relaxed ) )
        return;
    calls.fetch_add( 1, std::memory_order_relaxed );
    bytes.fetch_add( size, std::memory_order_relaxed );
}

/** @brief Reset and enable allocation recording. */
void begin() noexcept
{
    calls.store( 0, std::memory_order_relaxed );
    bytes.store( 0, std::memory_order_relaxed );
    enabled.store( true, std::memory_order_release );
}

/**
 * @brief Disable allocation recording and return the cumulative snapshot.
 * @return recorded C++ allocation calls and requested bytes
 */
AllocationSnapshot end() noexcept
{
    enabled.store( false, std::memory_order_release );
    return { calls.load( std::memory_order_relaxed ),
             bytes.load( std::memory_order_relaxed ) };
}

/**
 * @brief Allocate storage with the requested alignment.
 * @param size requested byte count
 * @param alignment requested alignment
 * @return allocated storage
 * @throws std::bad_alloc when allocation fails
 */
[[nodiscard]] void* allocate( std::size_t size, std::size_t alignment )
{
    record( size );
    void* storage = nullptr;
    if ( alignment <= alignof( std::max_align_t ) )
        storage = std::malloc( size == 0 ? 1 : size );
    else if ( ::posix_memalign( &storage, alignment, size == 0 ? 1 : size ) != 0 )
        storage = nullptr;
    if ( !storage )
        throw std::bad_alloc{};
    return storage;
}
} // namespace Feel::benchmarking::allocation_probe

/**
 * @brief Global scalar allocation hook used only by this benchmark executable.
 * @param size requested byte count
 * @return allocated storage
 */
void* operator new( std::size_t size )
{
    return Feel::benchmarking::allocation_probe::allocate(
        size, alignof( std::max_align_t ) );
}

/**
 * @brief Global array allocation hook used only by this benchmark executable.
 * @param size requested byte count
 * @return allocated storage
 */
void* operator new[]( std::size_t size )
{
    return Feel::benchmarking::allocation_probe::allocate(
        size, alignof( std::max_align_t ) );
}

/** @brief Release scalar storage allocated by the benchmark hook. */
void operator delete( void* storage ) noexcept
{
    std::free( storage );
}

/** @brief Release array storage allocated by the benchmark hook. */
void operator delete[]( void* storage ) noexcept
{
    std::free( storage );
}

/** @brief Release sized scalar storage allocated by the benchmark hook. */
void operator delete( void* storage, std::size_t ) noexcept
{
    std::free( storage );
}

/** @brief Release sized array storage allocated by the benchmark hook. */
void operator delete[]( void* storage, std::size_t ) noexcept
{
    std::free( storage );
}

/**
 * @brief Global aligned scalar allocation hook.
 * @param size requested byte count
 * @param alignment requested alignment
 * @return allocated aligned storage
 */
void* operator new( std::size_t size, std::align_val_t alignment )
{
    return Feel::benchmarking::allocation_probe::allocate(
        size, static_cast<std::size_t>( alignment ) );
}

/**
 * @brief Global aligned array allocation hook.
 * @param size requested byte count
 * @param alignment requested alignment
 * @return allocated aligned storage
 */
void* operator new[]( std::size_t size, std::align_val_t alignment )
{
    return Feel::benchmarking::allocation_probe::allocate(
        size, static_cast<std::size_t>( alignment ) );
}

/** @brief Release aligned scalar storage. */
void operator delete( void* storage, std::align_val_t ) noexcept
{
    std::free( storage );
}

/** @brief Release aligned array storage. */
void operator delete[]( void* storage, std::align_val_t ) noexcept
{
    std::free( storage );
}

/** @brief Release sized aligned scalar storage. */
void operator delete(
    void* storage, std::size_t, std::align_val_t ) noexcept
{
    std::free( storage );
}

/** @brief Release sized aligned array storage. */
void operator delete[](
    void* storage, std::size_t, std::align_val_t ) noexcept
{
    std::free( storage );
}

/**
 * @brief Nothrow scalar allocation hook.
 * @param size requested byte count
 * @return allocated storage or `nullptr`
 */
void* operator new( std::size_t size, std::nothrow_t const& ) noexcept
{
    try
    {
        return ::operator new( size );
    }
    catch ( ... )
    {
        return nullptr;
    }
}

/**
 * @brief Nothrow array allocation hook.
 * @param size requested byte count
 * @return allocated storage or `nullptr`
 */
void* operator new[]( std::size_t size, std::nothrow_t const& ) noexcept
{
    try
    {
        return ::operator new[]( size );
    }
    catch ( ... )
    {
        return nullptr;
    }
}

/** @brief Release nothrow scalar storage. */
void operator delete(
    void* storage, std::nothrow_t const& ) noexcept
{
    std::free( storage );
}

/** @brief Release nothrow array storage. */
void operator delete[](
    void* storage, std::nothrow_t const& ) noexcept
{
    std::free( storage );
}
