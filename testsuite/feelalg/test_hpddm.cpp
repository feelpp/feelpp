#define BOOST_TEST_MODULE test_hpddm

#include <feel/feelcore/testsuite.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feelalg/enums.hpp>
#include <feel/feelalg/preconditionerpetsc.hpp>
#include <feel/feelalg/topetsc.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelvf/vf.hpp>
#include <dlfcn.h>
#include <sstream>

using namespace Feel;
using namespace Feel::vf;

namespace
{

inline po::options_description makeOptions()
{
    po::options_description options( "test_hpddm options" );
    options.add( backend_options( "test-hpddm" ) );
    return options;
}

inline AboutData makeAbout()
{
    AboutData about( "test_hpddm",
                     "test_hpddm",
                     "0.1",
                     "test hpddm preconditioner plumbing",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2026 Feel++ Consortium" );

    about.addAuthor( "OpenAI", "developer", "support@openai.com", "" );
    return about;
}

std::string petscErrorSummary( int ierr )
{
    char const* text = nullptr;
    char const* specific = nullptr;
    PetscErrorMessage( ierr, &text, &specific );

    std::ostringstream os;
    os << "ierr=" << ierr
       << " text=" << ( text ? text : "<null>" )
       << " specific=" << ( specific ? specific : "<null>" );
    return os.str();
}

template <typename Signature>
Signature
lookupHpddmSymbol( char const* name )
{
    return reinterpret_cast<Signature>( dlsym( RTLD_DEFAULT, name ) );
}

} // namespace

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );

BOOST_AUTO_TEST_SUITE( hpddm_suite )

BOOST_AUTO_TEST_CASE( test_hpddm_string_conversion )
{
    BOOST_CHECK_EQUAL( pcTypeConvertStrToEnum( "hpddm" ), HPDDM_PRECOND );
}

BOOST_AUTO_TEST_CASE( test_hpddm_petsc_pc_type )
{
    if ( !lookupHpddmSymbol<void(*)()>( "PCCreate_HPDDM" ) )
    {
        BOOST_TEST_MESSAGE( "Skipping HPDDM PETSc runtime test: PETSc runtime does not export PCCreate_HPDDM" );
        BOOST_CHECK( true );
        return;
    }

    PC probe = nullptr;
    int ierr = PCCreate( Environment::worldComm().globalComm(), &probe );
    BOOST_REQUIRE_MESSAGE( ierr == 0 && probe != nullptr,
                           "PCCreate() failed while probing PETSc HPDDM runtime" );

    ierr = PetscPushErrorHandler( PetscReturnErrorHandler, nullptr );
    BOOST_REQUIRE_MESSAGE( ierr == 0,
                           "PetscPushErrorHandler() failed while probing PETSc HPDDM runtime" );

    ierr = PCSetType( probe, PCHPDDM );
    PetscPopErrorHandler();
    PETSc::PCDestroy( probe );

    BOOST_REQUIRE_MESSAGE( ierr == 0,
                           "PETSc exports PCCreate_HPDDM but PCSetType(pc, PCHPDDM) failed: "
                           + petscErrorSummary( ierr ) );

    auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );
    auto Xh = Pch<1>( mesh );
    auto b = backend();
    auto A = b->newMatrix( _trial=Xh, _test=Xh );
    auto u = Xh->element();
    auto v = Xh->element();
    form2( _trial=Xh, _test=Xh, _matrix=A ) =
        integrate( _range=elements( mesh ), _expr=idt( u )*id( v ) );
    A->close();

    auto Aaux = b->newMatrix( _trial=Xh, _test=Xh );
    form2( _trial=Xh, _test=Xh, _matrix=Aaux ) =
        integrate( _range=elements( mesh ), _expr=idt( u )*id( v ) );
    Aaux->close();

    auto p = Feel::preconditioner( _prefix="test-hpddm",
                                   _matrix=A,
                                   _pc=HPDDM_PRECOND,
                                   _backend=b,
                                   _worldcomm=Environment::worldCommPtr() );

    auto pp = toPETSc( p );
    BOOST_REQUIRE( pp );

    auto AauxPetsc = toPETSc( Aaux );
    BOOST_REQUIRE( AauxPetsc );

    PetscInt rowStart = 0, rowEnd = 0;
    ierr = MatGetOwnershipRange( AauxPetsc->mat(), &rowStart, &rowEnd );
    CHKERRABORT( Environment::worldComm().globalComm(), ierr );

    IS is = nullptr;
    ierr = ISCreateStride( Environment::worldComm().globalComm(), rowEnd-rowStart, rowStart, 1, &is );
    CHKERRABORT( Environment::worldComm().globalComm(), ierr );
    pp->attachHpddmAuxiliaryData( is, Aaux );
    ierr = ISDestroy( &is );
    CHKERRABORT( Environment::worldComm().globalComm(), ierr );

    using hpddm_set_auxiliary_mat_t = PetscErrorCode (*)( PC, IS, Mat, PetscErrorCode (*)(Mat, PetscReal, Vec, Vec, PetscReal, IS, void *), void * );
    BOOST_REQUIRE_MESSAGE( lookupHpddmSymbol<hpddm_set_auxiliary_mat_t>( "PCHPDDMSetAuxiliaryMat" ),
                           "PETSc exports PCCreate_HPDDM but does not export PCHPDDMSetAuxiliaryMat" );

    PetscInt auxiliaryIsRefBeforeInit = 0;
    ierr = PetscObjectGetReference( reinterpret_cast<PetscObject>( pp->hpddmAuxiliaryIS().get() ), &auxiliaryIsRefBeforeInit );
    CHKERRABORT( Environment::worldComm().globalComm(), ierr );

    PetscInt auxiliaryMatRefBeforeInit = 0;
    ierr = PetscObjectGetReference( reinterpret_cast<PetscObject>( AauxPetsc->mat() ), &auxiliaryMatRefBeforeInit );
    CHKERRABORT( Environment::worldComm().globalComm(), ierr );

    p->init();

    PC pc = pp->pc();
    BOOST_REQUIRE( pc );

    const char* pcType = nullptr;
    ierr = PCGetType( pc, &pcType );
    CHKERRABORT( Environment::worldComm().globalComm(), ierr );

    BOOST_REQUIRE( pcType != nullptr );
    BOOST_CHECK_EQUAL( std::string( pcType ), std::string( PCHPDDM ) );

    PetscInt auxiliaryIsRefAfterInit = 0;
    ierr = PetscObjectGetReference( reinterpret_cast<PetscObject>( pp->hpddmAuxiliaryIS().get() ), &auxiliaryIsRefAfterInit );
    CHKERRABORT( Environment::worldComm().globalComm(), ierr );
    BOOST_CHECK_GT( auxiliaryIsRefAfterInit, auxiliaryIsRefBeforeInit );

    PetscInt auxiliaryMatRefAfterInit = 0;
    ierr = PetscObjectGetReference( reinterpret_cast<PetscObject>( AauxPetsc->mat() ), &auxiliaryMatRefAfterInit );
    CHKERRABORT( Environment::worldComm().globalComm(), ierr );
    if ( Environment::worldComm().globalSize() == 1 )
        BOOST_CHECK_GT( auxiliaryMatRefAfterInit, auxiliaryMatRefBeforeInit );
    else
        BOOST_TEST_MESSAGE( "Skipping distributed auxiliary matrix refcount check: HPDDM stores a local extracted submatrix in parallel" );

    p->clear();

    PetscInt auxiliaryIsRefAfterClear = 0;
    ierr = PetscObjectGetReference( reinterpret_cast<PetscObject>( pp->hpddmAuxiliaryIS().get() ), &auxiliaryIsRefAfterClear );
    CHKERRABORT( Environment::worldComm().globalComm(), ierr );
    BOOST_CHECK_LT( auxiliaryIsRefAfterClear, auxiliaryIsRefAfterInit );

    p->init();
    pc = pp->pc();
    BOOST_REQUIRE( pc );

    ierr = PCGetType( pc, &pcType );
    CHKERRABORT( Environment::worldComm().globalComm(), ierr );
    BOOST_REQUIRE( pcType != nullptr );
    BOOST_CHECK_EQUAL( std::string( pcType ), std::string( PCHPDDM ) );

    PetscInt auxiliaryIsRefAfterReinit = 0;
    ierr = PetscObjectGetReference( reinterpret_cast<PetscObject>( pp->hpddmAuxiliaryIS().get() ), &auxiliaryIsRefAfterReinit );
    CHKERRABORT( Environment::worldComm().globalComm(), ierr );
    BOOST_CHECK_GT( auxiliaryIsRefAfterReinit, auxiliaryIsRefAfterClear );
}

BOOST_AUTO_TEST_CASE( test_hpddm_auxiliary_attachment_api )
{
    auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );
    auto Xh = Pch<1>( mesh );
    auto b = backend();
    auto A = b->newMatrix( _trial=Xh, _test=Xh );
    A->close();

    auto p = Feel::preconditioner( _prefix="test-hpddm-aux",
                                   _matrix=A,
                                   _pc=HPDDM_PRECOND,
                                   _backend=b,
                                   _worldcomm=Environment::worldCommPtr() );
    auto pp = toPETSc( p );
    BOOST_REQUIRE( pp );

    auto Aaux = b->newMatrix( _trial=Xh, _test=Xh );
    Aaux->close();

    PetscInt index = 0;
    IS is = nullptr;
    int ierr = ISCreateGeneral( Environment::worldComm().globalComm(), 1, &index, PETSC_COPY_VALUES, &is );
    CHKERRABORT( Environment::worldComm().globalComm(), ierr );

    pp->attachHpddmAuxiliaryData( is, Aaux );

    BOOST_CHECK( pp->hasHpddmAuxiliaryMatrix() );
    BOOST_CHECK( pp->hasHpddmAuxiliaryIS() );
    BOOST_CHECK_EQUAL( pp->hpddmAuxiliaryMatrix().get(), Aaux.get() );
    BOOST_CHECK_EQUAL( pp->hpddmAuxiliaryIS().get(), is );

    ierr = ISDestroy( &is );
    CHKERRABORT( Environment::worldComm().globalComm(), ierr );
    BOOST_CHECK( is == nullptr );

    PetscInt recoveredSize = 0;
    ierr = ISGetLocalSize( reinterpret_cast<IS>( pp->hpddmAuxiliaryIS().get() ), &recoveredSize );
    CHKERRABORT( Environment::worldComm().globalComm(), ierr );
    BOOST_CHECK_EQUAL( recoveredSize, 1 );

    PreconditionerPetsc<double> ppCopy( *pp );
    BOOST_CHECK( ppCopy.hasHpddmAuxiliaryMatrix() );
    BOOST_CHECK( ppCopy.hasHpddmAuxiliaryIS() );
    BOOST_CHECK_EQUAL( ppCopy.hpddmAuxiliaryMatrix().get(), Aaux.get() );
    BOOST_CHECK_EQUAL( ppCopy.hpddmAuxiliaryIS().get(), pp->hpddmAuxiliaryIS().get() );
}

BOOST_AUTO_TEST_SUITE_END()
