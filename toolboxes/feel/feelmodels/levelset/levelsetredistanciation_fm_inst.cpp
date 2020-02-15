
#include "levelsetbaseconfig.h"
#include <feel/feelmodels/levelset/levelsetbase.hpp>
#include <feel/feelmodels/levelset/levelsetredistanciation_fm.cpp>

namespace Feel {

template class LevelSetRedistanciationFM< 
    typename FeelModels::LEVELSETBASE_CLASS_TYPE::space_levelset_type
    >;

}
