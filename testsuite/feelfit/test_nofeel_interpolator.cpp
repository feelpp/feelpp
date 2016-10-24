#include <feel/feel.hpp>

#include <feel/feelfit/interpolator.hpp>



#if FEELPP_HAS_GSL
// to compare against GSL
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#endif

using namespace Feel;

int main(int argc, char **argv)
{
    Environment env( _argc=argc, _argv=argv,
                     _about=about(_name="test_nofeel_fit",
                                  _author="Feel++ Consortium",
                                  _email="feelpp-devel@feelpp.org"));
    
    std::ifstream infile("data.txt");
    double a, b;
    int i = 0;
    int size = 5;
    double gsl_x[size], gsl_y[size];
    for (int i = 0; i < 5; ++i )
    {
        infile >> a >> b;
        gsl_x[i] = a;
        gsl_y[i] = b;
    }
#if FEELPP_HAS_GSL
    // Cspline
    gsl_interp_accel *acc_spline = gsl_interp_accel_alloc ();
    //gsl_spline *spline = gsl_spline_alloc (gsl_interp_akima, size);
    gsl_spline *cspline = gsl_spline_alloc (gsl_interp_cspline, size);
    gsl_spline_init (cspline, gsl_x, gsl_y, size);
    // Akima
    gsl_interp_accel *acc_akima = gsl_interp_accel_alloc ();
    gsl_spline *akima = gsl_spline_alloc (gsl_interp_akima, size);
    gsl_spline_init (akima, gsl_x, gsl_y, size);
#endif
    auto po = Interpolator::New(P0,    "data.txt");
    auto p1 = Interpolator::New(P1,    "data.txt");
    auto cb = Interpolator::New(Spline,"data.txt");
    auto ak = Interpolator::New(Akima, "data.txt");
    std::ofstream ofile("data_interp.txt");
    ofile << "x\tP0\tP1\tSpline\tAkima\tGSL\n";
    double gsl_ak = 0., err_ak = 0.;
    double gsl_cs = 0., err_cs = 0.;
    for(double d = gsl_x[0]; d < gsl_x[i-1]; d+= 0.02)
    {
#if FEELPP_HAS_GSL
        gsl_cs = gsl_spline_eval (cspline, d, acc_spline);
        gsl_ak = gsl_spline_eval (akima, d, acc_akima);
        err_cs += std::fabs(gsl_cs-(*cb)(d));
        err_ak += std::fabs(gsl_ak-(*ak)(d));
#endif
        ofile << d << "\t" << (*po)(d) << "\t" << (*p1)(d) << "\t" << (*cb)(d) << "\t" << (*ak)(d) << "\t" 
#if FEELPP_HAS_GSL
              << gsl_cs << "\t"  
              << gsl_ak << "\t"  
#endif
              << std::endl;
    }
#if FEELPP_HAS_GSL
    gsl_spline_free (cspline);
    gsl_interp_accel_free (acc_spline);
    gsl_spline_free (akima);
    gsl_interp_accel_free (acc_akima);
    std::cout << "err_ak = " << err_ak << std::endl;
    std::cout << "err_cs = " << err_cs << std::endl;
#endif

    return 0;
}

