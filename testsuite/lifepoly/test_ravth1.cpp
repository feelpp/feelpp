#define BOOST_TEST_MODULE Raviar-Thomas polynomials test
// Boost.Test
#define USE_TEST 1
#include <boost/test/unit_test.hpp>
using boost::unit_test::test_suite;

#include <life/lifepoly/raviartthomas.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <life/lifepoly/Polynomial_set1.hpp>


#include <life/lifecore/life.hpp>
#include <life/lifemesh/simplex.hpp>
#include <life/lifemesh/simplexproduct.hpp>
#include <life/lifemesh/refentity.hpp>
#include <life/lifemesh/geond.hpp>




BOOST_AUTO_TEST_SUITE( ravth_testsuite )


BOOST_AUTO_TEST_CASE( ravth1 )
{
    using namespace Life;
    typedef RaviartThomas<0>::apply<3>::type rt0_type;
    typedef GeoND<3,Simplex<3, 1, 3> >::point_type point_type;

      GeoND<3,Simplex<3, 1, 3> > tetra;
    point_type V1; V1( 0 )=-1;V1( 1 )=-1;V1( 2 )=-1;
    point_type V2; V2( 0 )=1;V2( 1 )=-1;V2( 2 )=-1;
    point_type V3; V3( 0 )=-1;V3( 1 )=1;V3( 2 )=-1;
    point_type V4; V4( 0 )=-1;V4( 1 )=-1;V4( 2 )=1;
    tetra.setPoint( 0, V1 );
    tetra.setPoint( 1, V2 );
    tetra.setPoint( 2, V3 );
    tetra.setPoint( 3, V4 );
    tetra.update();

  std::cout << "[tetra] barycenter = " << tetra.barycenter() << "\n";
  std::cout<< "[tetra] measure = " << tetra.measure()<< std::endl ;
  std::cout<< tetra.faceMeasure(0) << std::endl ;
  std::cout<< tetra.faceMeasure(1) << std::endl ;
  std::cout<< tetra.faceMeasure(2) << std::endl ;
  std::cout<< tetra.faceMeasure(3) << std::endl ;




    rt0_type rt0;

  matrix_type m(12 ,4)   ;

  m(0 ,0)=m(0 ,1)=m(1 ,0)=m(1 ,2)=m(2 ,0)=m(2 ,3)=m(3 ,1)=m(4 ,0)
	 =m(4 ,2)=m(5 ,0)=m(5 ,3)=m(6,0) = m(6,1) =m(7 ,2) =m(8 , 0)
     =m(8, 3)=m(9,0) = m(9 ,1) =m(10,0) =m(10 ,2)
     =m(11,3) =  1./(3*tetra.measure()) ;
   m(0 ,2)=  m(0 ,3)= m(1 ,1)= m(1 ,3)= m(2 ,1)=m(2 ,2)
  =m(3 ,2)= m(3 ,3)= m(4 ,1)=  m(4 ,3)=m(5 ,1)= m(5 ,2)
        = m(6 ,2)= m(6 ,3) = m(7 ,1)= m(7 ,3)= m(8 ,1)= m(8 ,2)=
    m(9 ,2)= m(9 ,3)= m(10 ,1)= m(10 ,3)= m(11 ,1)= m(11 ,2)= 0. ;
   m(3 ,0)=m(7 ,0)=m(11 , 0) = -1./(3*tetra.measure()) ;

   Points_type p1(3); p1(0)=-1.; p1(1) =-1.; p1(2) = -1. ;
   /*  Points_type p2(3); p2(0)= 1.; p2(1) =-1.; p2(2) = -1. ;
   Points_type p3(3); p3(0)= -1.; p3(1) =1.; p3(2) = -1. ;
   /*  Points_type p4(3); p4(0)=-1.; p4(1) =-1.; p4(2) =  1. ;
   Points_type p5(3); p5(0)=-1.; p5(1) = 0.; p5(2) = -1. ;
   Points_type p6(3); p6(0)=-1.; p6(1) =-1.; p6(2) = 0.;
   Points_type p7(3); p7(0)= 0.; p7(1) = -1.; p7(2) = -1.;
   Points_type p8(3); p8(0)= -1./3.; p8(1) =-1./3.; p8(2) =-1./3.;
   Points_type p9(3); p9(0)=-1./3.; p9(1) = 1./3.;p9(2) = -1./3.;
   Points_type p10(3); p10(0)= 1./3.; p10(1) =-1./3.;p10(2) =-1./3.;*/
   // vectors_type P(10 ) ;
   vectors_type P(1 ) ;

   //   P(0)=p1;P(1)=p2;P(2)=p3;P(3)=p4;P(4)=p5;P(5)=p6;P(6)=p7;P(7)=p8;P(8)=p9;P(9)=p10 ;

         P(0)=p1;
   Polynomialset<20,12,4,1> Poly(m) ;

   std::cout << "**************************************************Analytic method*********************************************" << std::endl ;
   std::cout<<"coefficients= " << m << std::endl ;
   std::cout <<"Points= " << P << std::endl ;
   std::cout<<"evaluate_at_Points= " <<Poly.evaluates_Points ( P ) << std::endl ;


std::cout << "*****************************************Life method***********************************************************" << std::endl ;

       rt0_type::points_type pts(3,1);

   pts(0 ,0) = -1. ; pts(1 ,0) = -1. ; pts(2 ,0) = -1. ;
   /*  pts(0 ,0) =  1. ; pts(1 ,0) = -1. ; pts(2 ,0) = -1. ;
   pts(0 ,2) = -1. ; pts(1 ,2) =  1. ; pts(2 ,2) =-1. ;
   pts(0 ,3) = -1. ; pts(1 ,3) = -1. ; pts(2 ,3) = 1. ;
   pts(0 ,4) = -1.;  pts(1 ,4) = 0.  ; pts(2 ,4) = -1. ;
   pts(0 ,5) = -1.; pts(1 ,5) = -1.; pts(2 ,5) =  0. ;
   pts(0 ,6)= 0.; pts(1 ,6) =-1.;  pts(2 ,6) = -1. ;
   pts(0 ,7)=-1./3.; pts(1 ,7) = -1./3.;  pts(2 ,7) = -1./3. ;
   pts(0 ,8)= -1./3.; pts(1 ,8) =1./3.; pts(2 ,8) = -1./3. ;
   pts(0 ,9)= 1./3.; pts(1 ,9) =-1./3.; pts(2 ,9) = -1./3. ;
   */
/*
   Points_type p1(3); p1(0)=-1.; p1(1) =-1.; p1(2) = -1. ;
   Points_type p2(3); p2(0)= 1.; p2(1) =-1.; p1(2) = -1. ;
   Points_type p3(3); p3(0)= -1.; p3(1) =1.; p1(2) = -1. ;
   Points_type p4(3); p4(0)=-1.; p4(1) =-1.; p1(2) =  1. ;
   Points_type p5(3); p5(0)=-1.; p5(1) = 0.; p1(2) = -1. ;
   Points_type p6(3); p6(0)=-1.; p6(1) =-1.; p6(2) = 0.;
   Points_type p7(2); p7(0)= 0.; p7(1) = -1.; p7(2) = -1.;
   Points_type p8(2); p8(0)= -1./3.; p8(1) =-1./3.; p8(2) =-1./3.;
   Points_type p9(2); p9(0)=-1./3.; p9(1) = 1./3.;p9(2) = -1./3.;
   Points_type p10(2); p10(0)= 1./3.; p10(1) =-1./3.;p10(2) =-1./3.;
   vectors_type P(10 ) ;
*/



   std::cout << "pts= " << pts << "\n" ;
   auto eval_at_pts = rt0.evaluate( pts );
   std::cout << "eval_at_pts= "<< eval_at_pts << "\n" ;

   // BOOST_CHECK (ublas::norm_inf(Poly.evaluates_Points( P )-eval_at_pts) < 0.001 )  ;


}
BOOST_AUTO_TEST_SUITE_END()

