#define BOOST_TEST_MODULE Raviart-Thomas polynomials test
// Boost.Test
#define USE_TEST 1
#include <boost/test/unit_test.hpp>
using boost::unit_test::test_suite;

#include <life/lifepoly/raviartthomas.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <life/lifepoly/Polynomial_set.hpp>



BOOST_AUTO_TEST_SUITE( ravth2_testsuite )


BOOST_AUTO_TEST_CASE( ravth2 )
{
    Life::Assert::setLog( "test_ravth.log");
    using namespace Life;
    BOOST_TEST_MESSAGE( "start" );
    typedef RaviartThomas<1>::apply<2>::type rt0_type;

    rt0_type rt0;
    BOOST_TEST_MESSAGE( "RT_1 instantiated" );

    matrix_type m(16 ,6)   ;

    /* for(uint k = 0 ; k < 12 ;++ k )
       {
       m(k , 3) = m(k , 4) = m(k ,5) =0 ;
       }
    */
    m(0 ,0)=m(0 ,1)=m(1 ,0)=m(1 ,2)=m(2 ,1)=m(3 ,0)=m(3 ,2)
        =m(4 ,0)=m(4 ,1)=m(5 ,2)  = 1./4. ;
    m(2 ,0)=m(5 ,0)= -1./4. ;
    m(6 ,0)=m(6 ,1)= m(13 ,0)=m(14 ,0) = -3./4. ;
    m(7 ,0)=m(7 ,2)= 3./4. ;

    m(8 ,0)=m(8 ,1)= m(10 ,0)=m(10 ,1)= m(12 ,0)=m(15 ,0) = -6./4. ;

    m(9 ,0)=m(9 ,2)= m(11 ,0)=m(11 ,2)= 6./4. ;

    m(11 ,1)= 12./4. ;

    m(12 ,1)=m(15 ,2)= -8./4. ;
    m(12 ,2)=m(12 ,3)= m(15 ,3)= -4./4. ;
    m(12 ,4)=m(13 ,5)= m(14 ,4)=m(15 ,5)= -2./4. ;
    m(13 ,2) = -5./4. ;

    for(uint k = 0 ; k < 12 ;++ k )
    {
        m(k , 3) = m(k , 4) = m(k ,5) =0 ;
    }
    BOOST_TEST_MESSAGE( "coeff matrix created " );
    // std::cout<<" coefficients= " << m << std::endl ;

/*
  m(0 ,0)=m(0 ,1)=m(1 ,0)=m(1 ,2)=m(2 ,1)=m(3 ,0)=m(3 ,2)=m(4 ,0)
  =m(4 ,1)=m(5 ,2) = 1./4.;
  m(0 ,2)=m(1 ,1)=m(2 ,2)=m(3 ,1)=m(4 ,2)=m(5 ,1) = 0. ;
  m(2 ,0)=m(5 ,0)= -1./4. ;
*/

    Points_type p1(2); p1(0)=-1./3.; p1(1) =-1./3.;
    Points_type p2(2); p2(0)=-1.; p2(1) =-1.;
    Points_type p3(2); p3(0)= 1.; p3(1) =-1.;
    Points_type p4(2); p4(0)=-1.; p4(1) = 1.;
    Points_type p5(2); p5(0)=-1.; p5(1) = 0.;
    Points_type p6(2); p6(0)=-1./3.; p6(1) =-1./3.;
    Points_type p7(2); p7(0)= 0.; p7(1) = 0.;
    Points_type p8(2); p8(0)= 0.; p8(1) =-1.;
    Points_type p9(2); p9(0)=-1./2.; p9(1) = 1./2.;
    Points_type p10(2); p10(0)= 1./2.; p10(1) =-1./2.;
    vectors_type P(10 ) ;

    P(0)=p1;P(1)=p2;P(2)=p3;P(3)=p4;P(4)=p5;P(5)=p6;P(6)=p7;P(7)=p8;P(8)=p9;P(9)=p10 ;
    BOOST_TEST_MESSAGE( "points created " );

    Polynomialset<10,16,6,10> Poly(m) ;
    //  Points_type p(2); p(0)= 0.; p(1) = 0.;
    //  vectors_type P(1 ) ;
    //  P(0)=p ;
    BOOST_TEST_MESSAGE( "polyset created " );

    std::cout << "**************************************************Analytic method*********************************************" << std::endl ;
    std::cout<<"coefficients= " << m << std::endl ;
    std::cout <<"Points= " << P << std::endl ;
    std::cout<<"evaluate_at_Points= " <<Poly.evaluate_Points ( P ) << std::endl ;


    std::cout << "*****************************************Life method***********************************************************" << std::endl ;

    rt0_type::points_type pts(2,10);

    // pts(0 ,0) = 0. ; pts(1 ,0) = 0. ;

    pts(0 ,0) = -1./3. ; pts(1 ,0) = -1./3. ;
    pts(0 ,1) = -1. ; pts(1 ,1) = -1. ;
    pts(0 ,2) =  1. ; pts(1 ,2) = -1. ;
    pts(0 ,3) = -1. ; pts(1 ,3) =  1. ;
    pts(0 ,4) = -1. ; pts(1 ,4) =  0. ;
    pts(0 ,5)=-1./3.; pts(1 ,5) =-1./3.;
    pts(0 ,6)= 0.; pts(1 ,6) = 0.;
    pts(0 ,7)= 0.; pts(1 ,7) =-1.;
    pts(0 ,8)=-1./2.; pts(1 ,8) = 1./2.;
    pts(0 ,9)= 1./2.; pts(1 ,9) =-1./2.;

    std::cout << "pts= " << pts << "\n" ;
    auto eval_at_pts = rt0.evaluate( pts );
    std::cout << "eval_at_pts= "<< eval_at_pts << "\n" ;

    //  BOOST_CHECK (ublas::norm_inf(Poly.evaluate_Points( P )-eval_at_pts) < 1e-10 )  ;


}
BOOST_AUTO_TEST_SUITE_END()

