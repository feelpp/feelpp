/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*-

    This file is part of the Feel library

    Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
             Date: 2025-08-23

    Copyright (C) 2025 Université de Strasbourg
    SPDX-License-Identifier: LGPL-2.1-or-later
*/
#pragma once
#include <petsc.h>
#include <petscvec.h>
#include <stdexcept>

namespace Feel
{

inline MPI_Comm vecComm( Vec v )
{
    if ( !v ) 
        throw std::invalid_argument( "vecComm: null Vec pointer" );
    
    MPI_Comm comm = MPI_COMM_WORLD;
    PetscErrorCode ierr = PetscObjectGetComm( reinterpret_cast<PetscObject>( v ), &comm );
    if ( ierr != 0 )
        throw std::runtime_error( "vecComm: PetscObjectGetComm failed" );
    return comm;
}

/**
 * @brief RAII Guard class for reading PETSc vector arrays.
 * 
 * Automatically calls VecGetArrayRead in constructor and VecRestoreArrayRead 
 * in destructor to prevent resource leaks and crashes.
 * 
 * @code
 * Vec v;
 * PetscReadArrayGuard guard( v );
 * const PetscScalar* data = guard.data();
 * // Automatic cleanup when guard goes out of scope
 * @endcode 
 */
class PetscReadArrayGuard
{
public:
    [[nodiscard]] explicit PetscReadArrayGuard( Vec v )
        : M_vec( v )
    {
        if ( !v )
            throw std::invalid_argument( "PetscReadArrayGuard: null Vec pointer" );
            
        PetscErrorCode ierr = VecGetArrayRead( M_vec, &M_ptr );
        if ( ierr != 0 ) {
            CHKERRABORT( vecComm( M_vec ), ierr );
        }
        
        ierr = VecGetOwnershipRange( M_vec, &M_lo, &M_hi );
        if ( ierr != 0 ) {
            VecRestoreArrayRead( M_vec, &M_ptr ); // cleanup on error
            CHKERRABORT( vecComm( M_vec ), ierr );
        }
    }

    // Non-copyable
    PetscReadArrayGuard( const PetscReadArrayGuard& ) = delete;
    PetscReadArrayGuard& operator=( const PetscReadArrayGuard& ) = delete;

    // Moveable
    PetscReadArrayGuard( PetscReadArrayGuard&& other ) noexcept
        : M_vec( other.M_vec ), M_ptr( other.M_ptr ), M_lo( other.M_lo ), M_hi( other.M_hi )
    {
        other.M_vec = nullptr; 
        other.M_ptr = nullptr; 
        other.M_lo = 0; 
        other.M_hi = 0;
    }

    PetscReadArrayGuard& operator=( PetscReadArrayGuard&& other ) noexcept
    {
        if ( this != &other )
        {
            // Properly release current resources
            cleanup();
            
            // Move from other
            M_vec = other.M_vec; 
            M_ptr = other.M_ptr; 
            M_lo = other.M_lo; 
            M_hi = other.M_hi;
            
            // Reset other
            other.M_vec = nullptr; 
            other.M_ptr = nullptr; 
            other.M_lo = 0; 
            other.M_hi = 0;
        }
        return *this;
    }

    ~PetscReadArrayGuard() noexcept
    {
        cleanup();
    }

    [[nodiscard]] const PetscScalar* data() const noexcept { return M_ptr; }
    [[nodiscard]] PetscInt firstLocal() const noexcept { return M_lo; }
    [[nodiscard]] PetscInt lastLocal() const noexcept { return M_hi; } // exclusive
    [[nodiscard]] PetscInt localSize() const noexcept { return M_hi - M_lo; }
    [[nodiscard]] bool valid() const noexcept { return M_vec != nullptr && M_ptr != nullptr; }

private:
    void cleanup() noexcept
    {
        if ( M_vec && M_ptr )
        {
            PetscErrorCode ierr = VecRestoreArrayRead( M_vec, &M_ptr );
            if ( ierr != 0 ) {
                // In destructor, we can't throw, so use CHKERRABORT
                CHKERRABORT( vecComm( M_vec ), ierr );
            }
        }
        M_vec = nullptr;
        M_ptr = nullptr;
    }

    Vec M_vec = nullptr;
    const PetscScalar* M_ptr = nullptr;
    PetscInt M_lo = 0, M_hi = 0;
};

/**
 * @brief RAII Guard class for writing PETSc vector arrays.
 * 
 * @code
 * Vec v;
 * PetscWriteArrayGuard guard( v );
 * PetscScalar* data = guard.data();
 * // Modify data...
 * // Automatic cleanup when guard goes out of scope
 * @endcode
 */
class PetscWriteArrayGuard
{
public:
    [[nodiscard]] explicit PetscWriteArrayGuard( Vec v, bool auto_assembly = false )
        : M_vec( v ), M_auto_assembly( auto_assembly )
    {
        if ( !v )
            throw std::invalid_argument( "PetscWriteArrayGuard: null Vec pointer" );
            
        PetscErrorCode ierr = VecGetArray( M_vec, &M_ptr );
        if ( ierr != 0 ) {
            CHKERRABORT( vecComm( M_vec ), ierr );
        }
        
        ierr = VecGetOwnershipRange( M_vec, &M_lo, &M_hi );
        if ( ierr != 0 ) {
            VecRestoreArray( M_vec, &M_ptr ); // cleanup on error
            CHKERRABORT( vecComm( M_vec ), ierr );
        }
    }

    // Non-copyable
    PetscWriteArrayGuard( const PetscWriteArrayGuard& ) = delete;
    PetscWriteArrayGuard& operator=( const PetscWriteArrayGuard& ) = delete;

    // Moveable
    PetscWriteArrayGuard( PetscWriteArrayGuard&& other ) noexcept
        : M_vec( other.M_vec ), M_ptr( other.M_ptr ), M_lo( other.M_lo ), M_hi( other.M_hi ),
          M_auto_assembly( other.M_auto_assembly )
    {
        other.M_vec = nullptr; 
        other.M_ptr = nullptr; 
        other.M_lo = 0; 
        other.M_hi = 0;
        other.M_auto_assembly = false;
    }

    PetscWriteArrayGuard& operator=( PetscWriteArrayGuard&& other ) noexcept
    {
        if ( this != &other )
        {
            // Properly release current resources
            cleanup();
            
            // Move from other
            M_vec = other.M_vec; 
            M_ptr = other.M_ptr; 
            M_lo = other.M_lo; 
            M_hi = other.M_hi;
            M_auto_assembly = other.M_auto_assembly;
            
            // Reset other
            other.M_vec = nullptr; 
            other.M_ptr = nullptr; 
            other.M_lo = 0; 
            other.M_hi = 0;
            other.M_auto_assembly = false;
        }
        return *this;
    }

    ~PetscWriteArrayGuard() noexcept
    {
        cleanup();
    }

    [[nodiscard]] PetscScalar* data() const noexcept { return M_ptr; }
    [[nodiscard]] PetscInt firstLocal() const noexcept { return M_lo; }
    [[nodiscard]] PetscInt lastLocal() const noexcept { return M_hi; } // exclusive
    [[nodiscard]] PetscInt localSize() const noexcept { return M_hi - M_lo; }
    [[nodiscard]] bool valid() const noexcept { return M_vec != nullptr && M_ptr != nullptr; }

    // Manual assembly control
    void assembly()
    {
        if ( M_vec ) {
            PetscErrorCode ierr = VecAssemblyBegin( M_vec );
            CHKERRABORT( vecComm( M_vec ), ierr );
            ierr = VecAssemblyEnd( M_vec );
            CHKERRABORT( vecComm( M_vec ), ierr );
        }
    }

private:
    void cleanup() noexcept
    {
        if ( M_vec && M_ptr )
        {
            PetscErrorCode ierr = VecRestoreArray( M_vec, &M_ptr );
            if ( ierr != 0 ) {
                CHKERRABORT( vecComm( M_vec ), ierr );
            }
            
            if ( M_auto_assembly ) {
                ierr = VecAssemblyBegin( M_vec );
                if ( ierr == 0 ) {
                    VecAssemblyEnd( M_vec );
                }
            }
        }
        M_vec = nullptr;
        M_ptr = nullptr;
    }

    Vec M_vec = nullptr;
    PetscScalar* M_ptr = nullptr;
    PetscInt M_lo = 0, M_hi = 0;
    bool M_auto_assembly = false;
};

} // namespace Feel
