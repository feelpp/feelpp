#define USE_BOOST_TEST 1
#define BOOST_TEST_MODULE test_form_eigen
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/creategmshmesh.hpp>
#include <feel/feelfilters/domain.hpp>
#include <feel/feelvf/form.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelvf/vf.hpp>
#include <feel/feelfilters/exporter.hpp>
#include "feel/feelalg/matrixeigendense.hpp"

using namespace Feel;
using namespace Feel::vf;

#define T_CONV Hypercube
#define T_ORDER 1
#define T_DIM 2

inline
AboutData
makeAbout()
{
    AboutData about( "test_form_eigen" ,
                     "test_form_eigen" );
    return about;
}

FEELPP_ENVIRONMENT_WITH_OPTIONS( makeAbout(), feel_options() )
BOOST_AUTO_TEST_SUITE( eigen_suite )

    BOOST_AUTO_TEST_CASE( test_0 )
{
    using space_t = Pch_type<Mesh<T_CONV<T_DIM,T_ORDER,T_DIM>>,T_ORDER>;

    auto mesh = createGMSHMesh( _mesh=new Mesh<T_CONV<T_DIM,T_ORDER,T_DIM>>,
                                _desc=domain( _name=( boost::format( "%1%-%2%-%3%-%4%" )
                                                      % soption(_name="gmsh.domain.shape")
                                                      % T_DIM % T_ORDER
                                                      % T_CONV<T_DIM,T_ORDER,T_DIM>::type() ).str() ,
                                              _dim=T_DIM, _order=T_ORDER,
                                              _convex=T_CONV<T_DIM,T_ORDER,T_DIM>::type()  ) );

    auto backend_eigen = Backend<double>::build("eigen_dense", "",Environment::worldCommSeq() );

    for ( auto const& it : elements(mesh) )
    {
        auto submesh = createSubmesh( _mesh=mesh, _range=idedelements(mesh, it.id()), _worldcomm=Environment::worldCommSeq() );
        auto Ph = space_t::New( _mesh=submesh );

        auto w = Ph->element( "W" );
        auto fA = form2( _trial=Ph, _test=Ph, _backend=backend_eigen );


    }

}

BOOST_AUTO_TEST_SUITE_END()
