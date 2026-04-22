#define BOOST_TEST_MODULE test_hpddm_laplacian

#include <feel/feelcore/testsuite.hpp>

#include <feel/feelalg/backend.hpp>
#include <feel/feelalg/enums.hpp>
#include <feel/feelalg/preconditionerpetsc.hpp>
#include <feel/feelalg/topetsc.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelfilters/loadmesh.hpp>
#include <feel/feelvf/vf.hpp>
#include <dlfcn.h>

using namespace Feel;
using namespace Feel::vf;

namespace
{

constexpr char const* hpddmLaplacianPrefix = "test-hpddm-laplacian";

inline po::options_description makeOptions()
{
    po::options_description options( "test_hpddm_laplacian options" );
    options.add( backend_options( hpddmLaplacianPrefix ) );
    return options;
}

inline AboutData makeAbout()
{
    AboutData about( "test_hpddm_laplacian",
                     "test_hpddm_laplacian",
                     "0.1",
                     "test hpddm laplacian solve",
                     Feel::AboutData::License_GPL,
                     "Copyright (c) 2026 Feel++ Consortium" );

    about.addAuthor( "OpenAI", "developer", "support@openai.com", "" );
    return about;
}

template <typename Signature>
Signature
lookupHpddmSymbol( char const* name )
{
    return reinterpret_cast<Signature>( dlsym( RTLD_DEFAULT, name ) );
}

bool petscHasHpddmRuntime()
{
    if ( !lookupHpddmSymbol<void(*)()>( "PCCreate_HPDDM" ) ||
         !lookupHpddmSymbol<void(*)()>( "PCHPDDMSetAuxiliaryMat" ) )
        return false;

    PC pc = nullptr;
    int ierr = PCCreate( Environment::worldComm().globalComm(), &pc );
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
}

} // namespace

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), makeOptions() );

BOOST_AUTO_TEST_CASE( test_hpddm_laplacian_form_solve )
{
    if ( !petscHasHpddmRuntime() )
    {
        BOOST_TEST_MESSAGE( "Skipping HPDDM form solve smoke test: PETSc runtime does not provide HPDDM support" );
        BOOST_CHECK( true );
        return;
    }

    auto mesh = loadMesh( _mesh=new Mesh<Simplex<2>> );
    auto Vh = Pch<2>( mesh );
    auto u = Vh->element();
    auto v = Vh->element();

    auto exact = Px()*Px() + Py()*Py();

    auto l = form1( _test=Vh );
    l = integrate( _range=elements( mesh ), _expr=cst( -4.0 )*id( v ) );

    auto a = form2( _trial=Vh, _test=Vh );
    a = integrate( _range=elements( mesh ), _expr=inner( gradt( u ), grad( v ) ) );
    a.deferDirichlet();
    a += on( _range=boundaryfaces( mesh ), _rhs=l, _element=u, _expr=exact );

    a.solve( _rhs=l, _solution=u, _name=hpddmLaplacianPrefix, _rebuild=true );

    double l2 = normL2( _range=elements( mesh ), _expr=idv( u ) - exact );
    BOOST_TEST_MESSAGE( "HPDDM form solve L2 error = " + std::to_string( l2 ) );
    BOOST_CHECK_SMALL( l2, 1e-8 );
}
