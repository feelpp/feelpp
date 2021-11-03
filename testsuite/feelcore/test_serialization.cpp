#define USE_BOOST_TEST 1
#define BOOST_TEST_MODULE test_serialization

// need to be before call of BOOST_CLASS_EXPORT
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <feel/feelcore/testsuite.hpp>

#include <feel/feelcore/environment.hpp>
#include <feel/feelfilters/unitsquare.hpp>
#include <feel/feeldiscr/pch.hpp>
#include <feel/feelvf/vf.hpp>


using namespace Feel;
using namespace Feel::vf;

namespace test_serialization
{
void runTestVector( std::string const& format )
{
    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto v = Xh->element();

    auto v1 = backend()->newVector( Xh );
    auto v2 = backend()->newVector( Xh );

    auto ex = expr("sin(x):x");

    form1( _test=Xh, _vector=v1 )
        = integrate( _range=elements(mesh), _expr=ex*id(v) );
    v1->close();

    std::string file_name = "vector_archive";
    v1->save( file_name, format);
    v2->load( file_name, format );

    bool isOk=true;
    for (int k=0;k<Xh->nLocalDofWithGhost();++k )
        if ( std::abs((*v1)(k)-(*v2)(k))>1e-9 )
        {
            isOk = false;
            break;
        }
    BOOST_CHECK( isOk );

    v1->add(-1, v2);
    double err = v1->linftyNorm();
    BOOST_CHECK_SMALL(err, 1e-8);
}

}

FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( serialization_suite )

BOOST_AUTO_TEST_CASE( text_test )
{
    test_serialization::runTestVector( "text" );
}

BOOST_AUTO_TEST_CASE( xml_test )
{
    test_serialization::runTestVector( "xml" );
}

BOOST_AUTO_TEST_CASE( binary_test )
{
    test_serialization::runTestVector( "binary" );
}

BOOST_AUTO_TEST_CASE( binary_test_matrix )
{
    auto mesh = unitSquare();
    auto Xh = Pch<1>( mesh );
    auto v = Xh->element();

    auto m1 = backend()->newMatrix( _test=Xh, _trial=Xh );
    auto m2 = backend()->newMatrix( _test=Xh, _trial=Xh );

    form2( _test=Xh, _trial=Xh, _matrix=m1) =
        integrate( _range=elements(mesh), _expr=vf::sin(Px())*id(v)*idt(v) );
    m1->close();

    std::string fileBinary = "matrix_archive";
    m1->save( fileBinary );
    m2->load( fileBinary );

    m1->addMatrix(-1, m2);
    double err = m1->linftyNorm();
    BOOST_CHECK_SMALL(err, 1e-8);
}



BOOST_AUTO_TEST_SUITE_END()
