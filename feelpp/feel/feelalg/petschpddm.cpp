/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */

#include <feel/feelalg/petschpddm.hpp>
#include <dlfcn.h>

namespace Feel
{

void*
petscLookupHpddmSymbol( char const* name )
{
    return dlsym( RTLD_DEFAULT, name );
}

bool
petscHasHpddmSymbol( char const* name )
{
    return petscLookupHpddmSymbol( name ) != nullptr;
}

bool
petscHasHpddmAuxiliaryMatRuntime()
{
    return petscHasHpddmSymbol( "PCHPDDMSetAuxiliaryMat" );
}

bool
petscHasHpddmRuntime( MPI_Comm comm )
{
#if defined(PCHPDDM)
    if ( !petscHasHpddmSymbol( "PCCreate_HPDDM" ) ||
         !petscHasHpddmAuxiliaryMatRuntime() )
        return false;

    PC pc = nullptr;
    auto ierr = PCCreate( comm, &pc );
    if ( ierr != 0 || pc == nullptr )
        return false;

    ierr = PetscPushErrorHandler( PetscReturnErrorHandler, nullptr );
    if ( ierr != 0 )
    {
        PETSc::PCDestroy( pc );
        return false;
    }

    ierr = PCSetType( pc, PCHPDDM );
    PetscPopErrorHandler();
    PETSc::PCDestroy( pc );

    return ierr == 0;
#else
    return false;
#endif
}

} // namespace Feel
