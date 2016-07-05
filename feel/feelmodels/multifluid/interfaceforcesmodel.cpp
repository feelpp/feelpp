#include <feel/feelmodels/multifluid/interfaceforcesmodel.hpp>
#include <feel/feelmodels/multifluid/helfrichforcemodel.hpp>

template<typename LevelSetType>
InterfaceForcesModels<LevelSetType>::self_ptrtype 
InterfaceForcesModel<LevelSetType>::build( 
            std::string const& type, 
            std::string const& prefix, 
            levelset_ptrtype const& levelset )
{
    if( type == "helfrich" )
        return self_ptrtype( new HelfrichForceModel( prefix, levelset ) );
    else
    {
        CHECK(false) << "unknown interface forces model name (" << type << " requested)\n";
        return self_ptrtype();
    }
}
