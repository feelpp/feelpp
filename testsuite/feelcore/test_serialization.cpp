#define USE_BOOST_TEST 1
#define BOOST_TEST_MODULE test_serialization

// need to be before call of BOOST_CLASS_EXPORT
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <testsuite.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelvf/vf.hpp>


using namespace Feel;
using namespace Feel::vf;

FEELPP_ENVIRONMENT_NO_OPTIONS
BOOST_AUTO_TEST_SUITE( serialization_suite )

BOOST_AUTO_TEST_CASE( text_test )
{
    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto v = Xh->element();

    auto v1 = backend()->newVector( Xh );
    auto v2 = backend()->newVector( Xh );

    auto ex = expr("sin(x):x");

    form1( _test=Xh, _vector=v1) = integrate( elements(mesh), ex );

    std::string path= ".";
    v1->save( _path=path, _type="text");
    v2->load( _path=path, _type="text");

    v1->add(-1, v2);
    double err = v1->linftyNorm();
    BOOST_CHECK_SMALL(err, 1e-8);
}

BOOST_AUTO_TEST_CASE( xml_test )
{
    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto v = Xh->element();

    auto v1 = backend()->newVector( Xh );
    auto v2 = backend()->newVector( Xh );

    auto ex = expr("sin(x):x");

    form1( _test=Xh, _vector=v1) = integrate( elements(mesh), ex );

    std::string path= ".";
    v1->save( _path=path, _type="xml");
    v2->load( _path=path, _type="xml");

    v1->add(-1, v2);
    double err = v1->linftyNorm();
    BOOST_CHECK_SMALL(err, 1e-8);
}

BOOST_AUTO_TEST_CASE( binary_test )
{
    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto v = Xh->element();

    auto v1 = backend()->newVector( Xh );
    auto v2 = backend()->newVector( Xh );

    auto ex = expr("sin(x):x");

    form1( _test=Xh, _vector=v1) = integrate( elements(mesh), ex );

    std::string path= ".";
    v1->save( _path=path, _type="binary");
    v2->load( _path=path, _type="binary");

    v1->add(-1, v2);
    double err = v1->linftyNorm();
    BOOST_CHECK_SMALL(err, 1e-8);
}

#if defined(PETSC_HAVE_HDF5)
BOOST_AUTO_TEST_CASE( hdf5_test )
{
    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto v = Xh->element();

    auto v1 = backend()->newVector( Xh );
    auto v2 = backend()->newVector( Xh );

    auto ex = expr("sin(x):x");

    form1( _test=Xh, _vector=v1) = integrate( elements(mesh), ex );

    std::string path= ".";
    v1->save( _path=path, _type="hdf5");
    v2->load( _path=path, _type="hdf5");

    v1->add(-1, v2);
    double err = v1->linftyNorm();
    BOOST_CHECK_SMALL(err, 1e-8);
}
#endif
BOOST_AUTO_TEST_SUITE_END()
