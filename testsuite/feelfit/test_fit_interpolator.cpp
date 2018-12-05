
#define BOOST_TEST_MODULE test_fit_interpolator
#include <feel/feelcore/testsuite.hpp>

#include <feel/feelfit/interpolator.hpp>

#if defined( FEELPP_HAS_GSL )
// to compare against GSL
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#endif


FEELPP_ENVIRONMENT_NO_OPTIONS

BOOST_AUTO_TEST_SUITE( test_fit_interpolator_suite )

BOOST_AUTO_TEST_CASE( test_interpolator_p1 )
{
    using namespace Feel;
    std::vector<std::pair<double,double>> interpData;
    interpData.push_back( std::make_pair(0,0) );
    interpData.push_back( std::make_pair(1,2) );
    interpData.push_back( std::make_pair(2,1) );

    auto p1 = Interpolator::New(P1,interpData);
    for ( double d=0;d<=1;d=d+0.01 )
        BOOST_CHECK_SMALL( math::abs( 2.*d - p1->operator()(d) ), 1e-9);
    for ( double d=1;d<=2;d=d+0.01 )
        BOOST_CHECK_SMALL( math::abs( (3.-d) - p1->operator()(d) ), 1e-9);
}

#if defined( FEELPP_HAS_GSL )
BOOST_AUTO_TEST_CASE( test_interpolator_vs_gls )
{
    using namespace Feel;

    double a, b;
    int samplingSize = 50;//5;
    double gsl_x[samplingSize], gsl_y[samplingSize];
    std::vector<std::pair<double,double>> interpData;
    for ( int i = 0; i < samplingSize; ++i )
    {
        a=2.0*i;
        b=math::cos(3*a*a*a)+math::sqrt(2*a*a);
        interpData.push_back( std::make_pair(a,b) );
        gsl_x[i] = a;
        gsl_y[i] = b;
    }

    // linear
    gsl_interp_accel *acc_linear = gsl_interp_accel_alloc ();
    gsl_spline *linear = gsl_spline_alloc (gsl_interp_linear, samplingSize);
    gsl_spline_init( linear, gsl_x, gsl_y, samplingSize );
    // Cspline
    gsl_interp_accel *acc_spline = gsl_interp_accel_alloc();
    gsl_spline *cspline = gsl_spline_alloc( gsl_interp_cspline, samplingSize );
    gsl_spline_init( cspline, gsl_x, gsl_y, samplingSize );
    // Akima
    gsl_interp_accel *acc_akima = gsl_interp_accel_alloc();
    gsl_spline *akima = gsl_spline_alloc( gsl_interp_akima, samplingSize );
    gsl_spline_init( akima, gsl_x, gsl_y, samplingSize );

    auto po = Interpolator::New(P0,interpData);
    auto p1 = Interpolator::New(P1,interpData);
    auto cb = Interpolator::New(Spline,interpData);
    auto ak = Interpolator::New(Akima,interpData);

    for( double d = gsl_x[0]; d < gsl_x[samplingSize-1]; d+= 0.01 )
    {
        double gsl_p1 = gsl_spline_eval(linear, d, acc_linear);
        double gsl_cs = gsl_spline_eval(cspline, d, acc_spline);
        double gsl_ak = gsl_spline_eval(akima, d, acc_akima);
        double err_p1 = math::abs(gsl_p1-(*p1)(d));
        double err_cs = math::abs(gsl_cs-(*cb)(d));
        double err_ak = math::abs(gsl_ak-(*ak)(d));
        BOOST_CHECK_SMALL(err_p1, 1e-9);
        BOOST_CHECK_SMALL(err_cs, 1e-9);
        BOOST_CHECK_SMALL(err_ak, 1e-9);
    }

    gsl_spline_free(linear);
    gsl_interp_accel_free(acc_linear);
    gsl_spline_free(cspline);
    gsl_interp_accel_free(acc_spline);
    gsl_spline_free(akima);
    gsl_interp_accel_free(acc_akima);
}
#endif

BOOST_AUTO_TEST_SUITE_END()

