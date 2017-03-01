#define USE_BOOST_TEST 1
#define BOOST_TEST_MODULE test_serialization

// need to be before call of BOOST_CLASS_EXPORT
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <testsuite/testsuite.hpp>

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
    v1->close();


    std::string fileText = boost::str(boost::format("text_archive_%1%") % Environment::rank());

    {
        std::ofstream ofs(fileText);
        boost::archive::text_oarchive oa(ofs);
        oa << v1;
    }

    {
        std::ifstream ifs(fileText);
        boost::archive::text_iarchive ia(ifs);
        ia >> v2;
    }

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
    v1->close();

    std::string fileXml = boost::str(boost::format("xml_archive_%1%") % Environment::rank());

    {
        std::ofstream ofs(fileXml);
        boost::archive::xml_oarchive oa(ofs);
        oa << boost::serialization::make_nvp("vec", v1);
    }

    {
        std::ifstream ifs(fileXml);
        boost::archive::xml_iarchive ia(ifs);
        ia >> boost::serialization::make_nvp("vec", v2);
    }

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
    v1->close();


    std::string fileBinary = boost::str(boost::format("binary_archive_%1%") % Environment::rank());

    {
        std::ofstream ofs(fileBinary);
        boost::archive::binary_oarchive oa(ofs);
        oa << v1;
    }

    {
        std::ifstream ifs(fileBinary);
        boost::archive::binary_iarchive ia(ifs);
        ia >> v2;
    }

    v1->add(-1, v2);
    double err = v1->linftyNorm();
    BOOST_CHECK_SMALL(err, 1e-8);
}

BOOST_AUTO_TEST_CASE( binary_test_M )
{
    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto v = Xh->element();

    auto m1 = backend()->newMatrix( Xh, Xh );
    auto m2 = backend()->newMatrix( Xh, Xh );

    form2( _test=Xh, _trial=Xh, _matrix=m1) =
        integrate( elements(mesh), math::sin(Px())*id(v)*idt(v) );
    v1->close();


    std::string fileBinary = boost::str(boost::format("binary_archive_%1%") % Environment::rank());

    std::ofstream ofs(fileBinary);
    if (ofs)
    {
        boost::archive::binary_oarchive oa(ofs);
        oa << m1;
    }

    std::ifstream ifs(fileBinary);
    if (ifs)
    {
        boost::archive::binary_iarchive ia(ifs);
        ia >> m2;
    }

    m1->addMatrix(-1, m2);
    double err = m1->linftyNorm();
    BOOST_CHECK_SMALL(err, 1e-8);
}



BOOST_AUTO_TEST_SUITE_END()
