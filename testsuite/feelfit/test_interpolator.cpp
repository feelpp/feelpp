#include <feel/feelfit/interpolator.hpp>
#define GSL 1

#if GSL==1
// to compare against GSL
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#endif

int main(int argc, char **argv)
{
  std::vector<std::pair<double, double>> data;
  std::ifstream infile("data.txt");
  double a, b;
  int i = 0;
  int size = 5;
  double gsl_x[size], gsl_y[size];
  while (infile >> a >> b)
  {
    data.push_back({a,b});
    gsl_x[i] = a;
    gsl_y[i] = b;
    i++;
  }
#if GSL==1
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  //gsl_spline *spline = gsl_spline_alloc (gsl_interp_akima, size);
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, size);
  gsl_spline_init (spline, gsl_x, gsl_y, size);
#endif
  Interpol* po = Interpol::New(P0,data);
  Interpol* p1 = Interpol::New(P1,data);
  Interpol* cb = Interpol::New(Spline,data);
  Interpol* ak = Interpol::New(Akima,data);
  std::ofstream ofile("data_interp.txt");
  ofile << "x\tP0\tP1\tSpline\tAkima\tGSL\n";
  for(double d = 0.3; d < 6; d+= 0.02)
  {
#if GLS==1
    b = gsl_spline_eval (spline, d, acc);
#endif
    ofile << d << "\t" << (*po)(d) << "\t" << (*p1)(d) << "\t" << (*cb)(d) << "\t" << (*ak)(d) << "\t" 
#if GLS==1
      << b 
#endif
      << std::endl;
  }
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  return 0;
}

