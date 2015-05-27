

#include <feel/feelmodels2/thermodyn/thermodynamics.cpp>

namespace Feel {
namespace FeelModels {

#include "thermodynconfig.h"
template class ThermoDynamics< Simplex<THERMODYNAMICS_DIM,THERMODYNAMICS_ORDERGEO,THERMODYNAMICS_DIM>, THERMODYNAMICS_ORDERPOLY >;

}
}
