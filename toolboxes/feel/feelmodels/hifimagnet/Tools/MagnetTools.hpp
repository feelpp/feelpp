#ifndef __MagnetTools_Interface_H_
#define __MagnetTools_Interface_H_

#include <string>

// **** Add specific headers and defs for MagnetTools/B_Map *** //
bool omp_on = false;
bool omp_already_on = false;

int gsl_error_mode = 0;
int plot_mode = 0;
int check_mode = 0;
int sym_mode = 0;
int read_mode = 0;
int verbose_mode = 0;
int geom_mode = 0;
int ensight_on = 0;

std::string ensight_mode = "C Binary";

double T_ref = 20;
double factor_J = 1.e-08;
double factor_P = 1.e-06;
double factor_Sigma = 1.e-08;
double factor_Constraint = 1.e-08;
double units_by_defaults = 1;

#include "Constraints.h"
double ep = 0;
Constraints Actual_Constraints;

#include "Aubert_files.h"
#include "parallel.h"
#include "NumIntegration.h"
#include "eps_params.h"
int num_integ = 50;
int num_eval = 100;

// To avoid loading/writing an eps_params.dat
// this should be changed to allow use of bmap_options (see eps_params.h)
epsparams Precision(1.e-12, 1.e14, 1.e12, 1.e-14, 1.e-15);
gsl_integration_data ** Integ_ptr = NULL;

//#include <gsl/gsl_cadna.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_trig.h>

#endif
