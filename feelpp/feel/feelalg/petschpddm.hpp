/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

#ifndef FEELPP_ALG_PETSCHPDDM_HPP
#define FEELPP_ALG_PETSCHPDDM_HPP 1

#include <feel/feelcore/feelpetsc.hpp>

namespace Feel
{

FEELPP_EXPORT void* petscLookupHpddmSymbol( char const* name );

template <typename Signature>
Signature
petscHpddmSymbol( char const* name )
{
    return reinterpret_cast<Signature>( petscLookupHpddmSymbol( name ) );
}

FEELPP_EXPORT bool petscHasHpddmSymbol( char const* name );
FEELPP_EXPORT bool petscHasHpddmAuxiliaryMatRuntime();
FEELPP_EXPORT bool petscHasHpddmRuntime( MPI_Comm comm );

} // namespace Feel

#endif
