

#include <feel/feelmodels2/thermodyn/thermodynbase.cpp>

namespace Feel {
namespace FeelModels {

#include "thermodynconfig.h"
template class ThermoDynamicsBase< Simplex<THERMODYNAMICS_DIM,THERMODYNAMICS_ORDERGEO,THERMODYNAMICS_DIM>, THERMODYNAMICS_ORDERPOLY >;

}
}
