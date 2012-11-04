
#define BOOST_TEST_MODULE kdtree tests
#include <testsuite/testsuite.hpp>


#include <feel/feelmesh/kdtree.hpp>
using namespace Feel;

FEELPP_ENVIRONMENT_NO_OPTIONS

bool
check_point( const KDTree::node_type & p1,const KDTree::node_type & p2 )
{
    bool b=true;

    for ( size_type i=0; i<p1.size(); ++i )
        if ( ( ( p1( i )-p2( i ) )*( p1( i )-p2( i ) ) ) > 1e-15 ) b=false;

    return b;
}

bool
check_inlist( const KDTree::points_search_type & L, const KDTree::node_type & p )
{
    KDTree::points_search_const_iterator it = L.begin();
    KDTree::points_search_const_iterator it_end = L.end();
    bool find=false;

    while ( it!=it_end && !find )
    {
        if ( check_point( boost::get<0>( *it ),p ) ) find=true;

        else ++it;
    }

    return find;
}

BOOST_AUTO_TEST_SUITE( kdtree_testsuite )

BOOST_AUTO_TEST_CASE( test_1d_kdtree )
{
    KDTree kd_tree;
    KDTree::node_type n( 1 );

    for ( double i=0; i<=2; i=i+0.01 )
    {
        n( 0 )=i;
        kd_tree.addPoint( n );
    }

    /*********************
     * case 1 :
     *********************/

    //the point of research
    KDTree::node_type p1( 1 );
    p1( 0 )=0.314;
    //the solution
    KDTree::node_type p11( 1 );
    p11( 0 )=0.31;
    KDTree::node_type p12( 1 );
    p12( 0 )=0.32;

    kd_tree.nbNearNeighbor( 2 );
    kd_tree.search( p1 );
    //  kd_tree.showResultSearch();
    KDTree::points_search_type L1=kd_tree.pointsNearNeighbor();

    BOOST_CHECK( check_point( boost::get<0>( L1[0] ),p11 ) );
    BOOST_CHECK( check_point( boost::get<0>( L1[1] ),p12 ) );

    /*********************
     * case 2 :
     *********************/
    KDTree::node_type p2( 1 );
    p2( 0 )=1.81;

    KDTree::node_type p21( 1 );
    p21( 0 )=1.80;
    KDTree::node_type p22( 1 );
    p22( 0 )=1.82;



    kd_tree.nbNearNeighbor( 3 );
    kd_tree.search( p2 );
    KDTree::points_search_type L2=kd_tree.pointsNearNeighbor();
    //kd_tree.showResultSearch();

    //the first is the point p2
    BOOST_CHECK( boost::get<4>( L2[0] )<1e-15 );
    //the others are all at equal distances
    BOOST_CHECK( math::pow( boost::get<4>( L2[1] )-boost::get<4>( L2[2] ),2 )<1e-15 );
    //verify that it's the good points
    BOOST_CHECK( check_point( boost::get<0>( L2[0] ),p2 ) );
    BOOST_CHECK( check_inlist( L2,p21 ) );
    BOOST_CHECK( check_inlist( L2,p22 ) );

}

BOOST_AUTO_TEST_CASE( test_2d_kdtree )
{
    KDTree kd_tree;
    KDTree::node_type n( 2 );

    for ( double i=0; i<=2; i=i+0.1 )
        for ( double j=0; j<=2; j=j+0.1 )
        {
            n( 0 )=i;
            n( 1 )=j;
            kd_tree.addPoint( n );
        }

    /*********************
     * case 1 :
     *********************/

    //the point of research
    KDTree::node_type p1( 2 );
    p1( 0 )=0.43;
    p1( 1 )=1.32;
    //the solution
    KDTree::node_type p11( 2 );
    p11( 0 )=0.4;
    p11( 1 )=1.3;
    KDTree::node_type p12( 2 );
    p12( 0 )=0.5;
    p12( 1 )=1.3;
    KDTree::node_type p13( 2 );
    p13( 0 )=0.4;
    p13( 1 )=1.4;
    KDTree::node_type p14( 2 );
    p14( 0 )=0.5;
    p14( 1 )=1.4;

    kd_tree.nbNearNeighbor( 4 );
    kd_tree.search( p1 );
    //  kd_tree.showResultSearch();
    KDTree::points_search_type L1=kd_tree.pointsNearNeighbor();

    BOOST_CHECK( check_point( boost::get<0>( L1[0] ),p11 ) );
    BOOST_CHECK( check_point( boost::get<0>( L1[1] ),p12 ) );
    BOOST_CHECK( check_point( boost::get<0>( L1[2] ),p13 ) );
    BOOST_CHECK( check_point( boost::get<0>( L1[3] ),p14 ) );

    /*********************
     * case 2 :
     *********************/
    KDTree::node_type p2( 2 );
    p2( 0 )=0.8;
    p2( 1 )=0.2;

    KDTree::node_type p21( 2 );
    p21( 0 )=0.8;
    p21( 1 )=0.1;
    KDTree::node_type p22( 2 );
    p22( 0 )=0.8;
    p22( 1 )=0.3;
    KDTree::node_type p23( 2 );
    p23( 0 )=0.7;
    p23( 1 )=0.2;
    KDTree::node_type p24( 2 );
    p24( 0 )=0.9;
    p24( 1 )=0.2;

    kd_tree.nbNearNeighbor( 5 );
    kd_tree.search( p2 );
    KDTree::points_search_type L2=kd_tree.pointsNearNeighbor();
    //kd_tree.showResultSearch();

    //the first is the point p2
    BOOST_CHECK( boost::get<4>( L2[0] )<1e-15 );
    //the others are all at equal distances
    BOOST_CHECK( math::pow( boost::get<4>( L2[1] )-boost::get<4>( L2[2] ),2 )<1e-15 );
    BOOST_CHECK( math::pow( boost::get<4>( L2[1] )-boost::get<4>( L2[3] ),2 )<1e-15 );
    BOOST_CHECK( math::pow( boost::get<4>( L2[1] )-boost::get<4>( L2[4] ),2 )<1e-15 );
    //verify that it's the good points
    BOOST_CHECK( check_point( boost::get<0>( L2[0] ),p2 ) );
    BOOST_CHECK( check_inlist( L2,p21 ) );
    BOOST_CHECK( check_inlist( L2,p22 ) );
    BOOST_CHECK( check_inlist( L2,p23 ) );
    BOOST_CHECK( check_inlist( L2,p24 ) );

}

BOOST_AUTO_TEST_CASE( test_3d_kdtree )
{

    KDTree kd_tree;
    KDTree::node_type n( 3 );

    for ( double i=0; i<=2; i=i+0.1 )
        for ( double j=0; j<=2; j=j+0.1 )
            for ( double k=0; k<=2; k=k+0.1 )
            {
                n( 0 )=i;
                n( 1 )=j;
                n( 2 )=k;
                kd_tree.addPoint( n );
            }

    /*********************
     * case 1 :
     *********************/

    //the point of research
    KDTree::node_type p1( 3 );
    p1( 0 )=0.43;
    p1( 1 )=1.32;
    p1( 2 )=1.835;

    //the solution
    KDTree::node_type p11( 3 );
    p11( 0 )=0.4;
    p11( 1 )=1.3;
    p11( 2 )=1.8;
    KDTree::node_type p12( 3 );
    p12( 0 )=0.5;
    p12( 1 )=1.3;
    p12( 2 )=1.8;
    KDTree::node_type p13( 3 );
    p13( 0 )=0.4;
    p13( 1 )=1.4;
    p13( 2 )=1.8;
    KDTree::node_type p14( 3 );
    p14( 0 )=0.5;
    p14( 1 )=1.4;
    p14( 2 )=1.8;
    KDTree::node_type p15( 3 );
    p15( 0 )=0.4;
    p15( 1 )=1.3;
    p15( 2 )=1.9;
    KDTree::node_type p16( 3 );
    p16( 0 )=0.5;
    p16( 1 )=1.3;
    p16( 2 )=1.9;
    KDTree::node_type p17( 3 );
    p17( 0 )=0.4;
    p17( 1 )=1.4;
    p17( 2 )=1.9;
    KDTree::node_type p18( 3 );
    p18( 0 )=0.5;
    p18( 1 )=1.4;
    p18( 2 )=1.9;

    kd_tree.nbNearNeighbor( 8 );
    kd_tree.search( p1 );
    //kd_tree.showResultSearch();
    KDTree::points_search_type L1=kd_tree.pointsNearNeighbor();

    //on assure la liste croissante

    double var_temp=0;

    for ( size_type i=0; i<8; ++i )
    {
        BOOST_CHECK( var_temp <= boost::get<4>( L1[i] ) );
        var_temp = boost::get<4>( L1[i] );
    }

    //on s'assure que les point y sont
    BOOST_CHECK( check_inlist( L1,p11 ) );
    BOOST_CHECK( check_inlist( L1,p12 ) );
    BOOST_CHECK( check_inlist( L1,p13 ) );
    BOOST_CHECK( check_inlist( L1,p14 ) );
    BOOST_CHECK( check_inlist( L1,p15 ) );
    BOOST_CHECK( check_inlist( L1,p16 ) );
    BOOST_CHECK( check_inlist( L1,p17 ) );
    BOOST_CHECK( check_inlist( L1,p18 ) );


    /*********************
     * case 2 :
     *********************/

    KDTree::node_type p2( 3 );
    p2( 0 )=0.9;
    p2( 1 )=1.4;
    p2( 2 )=1.0;
    //the solution
    KDTree::node_type p21( 3 );
    p21( 0 )=0.9;
    p21( 1 )=1.3;
    p21( 2 )=1.0;
    KDTree::node_type p22( 3 );
    p22( 0 )=0.9;
    p22( 1 )=1.5;
    p22( 2 )=1.0;
    KDTree::node_type p23( 3 );
    p23( 0 )=0.8;
    p23( 1 )=1.4;
    p23( 2 )=1.0;
    KDTree::node_type p24( 3 );
    p24( 0 )=1.0;
    p24( 1 )=1.4;
    p24( 2 )=1.0;
    KDTree::node_type p25( 3 );
    p25( 0 )=0.9;
    p25( 1 )=1.4;
    p25( 2 )=1.1;
    KDTree::node_type p26( 3 );
    p26( 0 )=0.9;
    p26( 1 )=1.4;
    p26( 2 )=0.9;

    kd_tree.nbNearNeighbor( 7 );
    kd_tree.search( p2 );
    KDTree::points_search_type L2=kd_tree.pointsNearNeighbor();
    // kd_tree.showResultSearch();

    //the first is the point p2
    BOOST_CHECK( boost::get<4>( L2[0] )<1e-15 );
    //the others are all at equal distances
    BOOST_CHECK( math::pow( boost::get<4>( L2[1] )-boost::get<4>( L2[2] ),2 )<1e-15 );
    BOOST_CHECK( math::pow( boost::get<4>( L2[1] )-boost::get<4>( L2[3] ),2 )<1e-15 );
    BOOST_CHECK( math::pow( boost::get<4>( L2[1] )-boost::get<4>( L2[4] ),2 )<1e-15 );
    BOOST_CHECK( math::pow( boost::get<4>( L2[1] )-boost::get<4>( L2[5] ),2 )<1e-15 );
    BOOST_CHECK( math::pow( boost::get<4>( L2[1] )-boost::get<4>( L2[6] ),2 )<1e-15 );

    //verify that it's the good points
    BOOST_CHECK( check_point( boost::get<0>( L2[0] ),p2 ) );
    BOOST_CHECK( check_inlist( L2,p21 ) );
    BOOST_CHECK( check_inlist( L2,p22 ) );
    BOOST_CHECK( check_inlist( L2,p23 ) );
    BOOST_CHECK( check_inlist( L2,p24 ) );
    BOOST_CHECK( check_inlist( L2,p25 ) );
    BOOST_CHECK( check_inlist( L2,p26 ) );


}


BOOST_AUTO_TEST_SUITE_END() // kdtree test suite


