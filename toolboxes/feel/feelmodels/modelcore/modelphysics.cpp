/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4
 */

#include <feel/feelmodels/modelcore/modelphysics.hpp>

namespace Feel
{
namespace FeelModels
{

template <uint16_type Dim>
ModelPhysic<Dim>::ModelPhysic( std::string const& type, std::string const& name, ModelModel const& model )
    :
    M_type( type ),
    M_name( name )
{

    // save name of submodel but not yet intialized
    for ( std::string const& submodel : model.submodels() )
        M_subphysics[submodel];

    if ( type == "GenericPDE" )
        return;

    if ( M_type == "thermo-electric" )
        M_subphysicsTypes.insert( { "heat","electric" } );
    else if ( M_type == "heat-fluid" )
        M_subphysicsTypes.insert( { "heat","fluid" } );


    material_property_shape_dim_type scalarShape = std::make_pair(1,1);
    material_property_shape_dim_type matrixShape = std::make_pair(nDim,nDim);

    this->addMaterialPropertyDescription( "density", "rho", { scalarShape } );
    if ( M_type == "heat" ||  M_type == "thermo-electric" || M_type == "heat-fluid" )
    {
        this->addMaterialPropertyDescription( "specific-heat-capacity", "Cp", { scalarShape } );
        this->addMaterialPropertyDescription( "thermal-expansion", "beta", { scalarShape } );
        this->addMaterialPropertyDescription( "thermal-conductivity", "k", { scalarShape,matrixShape } );
    }
    if ( M_type == "electric" || M_type == "thermo-electric" )
    {
        this->addMaterialPropertyDescription( "electric-conductivity", "sigma", { scalarShape } );
    }
}

template <uint16_type Dim>
std::shared_ptr<ModelPhysic<Dim>>
ModelPhysic<Dim>::New( std::string const& type, std::string const& name, ModelModel const& model )
{
    if ( type == "fluid" )
        return std::make_shared<ModelPhysicFluid<Dim>>( name, model );
    else
        return std::make_shared<ModelPhysic<Dim>>( type, name, model );
}

template <uint16_type Dim>
ModelPhysicFluid<Dim>::ModelPhysicFluid( std::string const& name, ModelModel const& model )
    :
    super_type( "fluid", name, model ),
    M_equation( "Navier-Stokes" )
{
    auto const& pt = model.ptree();
    if ( auto eq = pt.template get_optional<std::string>("equations") )
        this->setEquation( *eq );
    if ( auto eq = pt.template get_optional<std::string>("equation") )
        this->setEquation( *eq );
}

template <uint16_type Dim>
void
ModelPhysicFluid<Dim>::setEquation( std::string const& eq )
{
    CHECK( eq == "Navier-Stokes" || eq == "Stokes" || eq == "StokesTransient" ) << "invalid equation of fluid : " << eq;
    M_equation = eq;
}

template <uint16_type Dim>
void
ModelPhysics<Dim>::initPhysics( std::string const& name, ModelModels const& models, subphysic_description_type const& subPhyicsDesc )
{
    std::string const& type = M_physicType;

    M_physicDefault = name;
    auto const& theGlobalModel = models.model( name );
    M_physics.emplace( name, ModelPhysic<Dim>::New( type, name, theGlobalModel ) );
    for ( auto const& [variantName,variantModel] : theGlobalModel.variants() )
        M_physics.emplace( variantName, ModelPhysic<Dim>::New( type, variantName, variantModel ) );

    // list of subphysics from description
    std::map<std::string, std::set<std::string>> mapSubPhysicsTypeToDefaultNames;
    std::vector<std::shared_ptr<ModelPhysic<Dim>>> parentPhysicToUpdate;
    for ( auto const& [physicName,physicData] : M_physics )
    {
        for ( std::string const& subPhysicType : physicData->subphysicsTypes() )
        {
            auto itFindSubphysic = std::find_if( physicData->subphysics().begin(), physicData->subphysics().end(), [&subPhysicType]( auto const& e ) {
                    if ( !e.second )
                        return false;
                    return e.second->type() == subPhysicType;
                } );
            if ( itFindSubphysic != physicData->subphysics().end() )
                continue;

            auto itFindTypeInDesc = subPhyicsDesc.find( subPhysicType );
            CHECK( itFindTypeInDesc != subPhyicsDesc.end() ) << "type not given in subPhyicsDesc";
            std::string const& subPhysicDefaultName = itFindTypeInDesc->second;
            mapSubPhysicsTypeToDefaultNames[subPhysicType].insert( subPhysicDefaultName );
            parentPhysicToUpdate.push_back( physicData );
        }
    }

    // create subphysics
    bool useModelNameInJson = models.useModelName();
    for ( auto const& [subPhysicType,subPhysicDefaultNames] : mapSubPhysicsTypeToDefaultNames )
    {
        for ( std::string const& subName : subPhysicDefaultNames )
        {
            auto const& theSubModel = useModelNameInJson && models.hasModel( subName )? models.model( subName ) : ModelModel{};
            M_physics.emplace( subName, ModelPhysic<Dim>::New( subPhysicType, subName, theSubModel ) );
            for ( auto const& [variantName,variantModel] : theSubModel.variants() )
                M_physics.emplace( variantName, ModelPhysic<Dim>::New( subPhysicType, variantName, variantModel ) );
        }
    }

    // update subphysics connection
    for ( auto const& parentPhysic : parentPhysicToUpdate )
    {
        std::string const& parentType = parentPhysic->type();
        std::string const& parentName = parentPhysic->name();

        // if subphysics defined but not init
        auto const& subphysicRegistered = parentPhysic->subphysics();
        for ( auto const& subPhysicPair : subphysicRegistered )
        {
            std::string subPhysicName = subPhysicPair.first;
            auto const& subPhysicData = subPhysicPair.second;
            if ( subPhysicData )
                continue;
            auto itFindSubphysic = std::find_if( M_physics.begin(), M_physics.end(), [&subPhysicName]( auto const& e ) { return /*e.second->type() == subPhysicType &&*/ e.second->name() == subPhysicName; } );
            CHECK( itFindSubphysic != M_physics.end() ) << "subphysic not found";
            parentPhysic->addSubphysic( itFindSubphysic->second );
        }

        // init all subphysics
        for ( std::string const& subPhysicType : parentPhysic->subphysicsTypes() )
        {
            auto itFindType = std::find_if( subphysicRegistered.begin(), subphysicRegistered.end(), [&subPhysicType]( auto const& e ) { return e.second->type() == subPhysicType; } );
            if ( itFindType != subphysicRegistered.end() )
                continue;

            auto itFindTypeInDesc = subPhyicsDesc.find( subPhysicType );
            CHECK( itFindTypeInDesc != subPhyicsDesc.end() ) << "type not given in subPhyicsDesc";
            std::string const& subPhysicDefaultName = itFindTypeInDesc->second;

            auto itFindSubphysic = std::find_if( M_physics.begin(), M_physics.end(), [&subPhysicType,&subPhysicDefaultName]( auto const& e ) { return e.second->type() == subPhysicType && e.second->name() == subPhysicDefaultName; } );
            CHECK( itFindSubphysic != M_physics.end() ) << "subphysic not found";
            parentPhysic->addSubphysic( itFindSubphysic->second );

        }
    }

}


template <uint16_type Dim>
std::set<std::string>
ModelPhysics<Dim>::physicsAvailable() const
{
    std::set<std::string> res;
    for ( auto const& [physicName,physicData] : M_physics )
        res.insert( physicName );
    return res;
}

template <uint16_type Dim>
std::set<std::string>
ModelPhysics<Dim>::physicsAvailable( std::string const& type ) const
{
    std::set<std::string> res;
    for ( auto const& [physicName,physicData] : M_physics )
        if ( physicData->type() == type )
            res.insert( physicName );
    return res;
}

template <uint16_type Dim>
std::set<std::string>
ModelPhysics<Dim>::physicsShared( std::string const& pname ) const
{
    std::set<std::string> res;
    auto itFindPhysic = M_physics.find( pname );
    if ( itFindPhysic != M_physics.end() )
        itFindPhysic->second->physicsShared( res );
    return res;
}

template <uint16_type Dim>
std::map<std::string,std::shared_ptr<ModelPhysic<Dim>>>
    ModelPhysics<Dim>::physics( std::string const& type ) const
{
    std::map<std::string,std::shared_ptr<ModelPhysic<nDim>>> res;
    for ( auto const& [physicName,physicData] : M_physics )
    {
        if ( physicData->type() == type )
            res[physicName] = physicData;
    }
    return res;
}

template <uint16_type Dim>
void
ModelPhysics<Dim>::setPhysics( std::map<std::string,std::shared_ptr<ModelPhysic<Dim>>> const& thePhysics, std::string const& physicDefault )
{
    M_physics.insert( thePhysics.begin(), thePhysics.end() );

    if ( !physicDefault.empty() )
        M_physicDefault = physicDefault;
}

template class ModelPhysic<2>;
template class ModelPhysic<3>;
template class ModelPhysics<2>;
template class ModelPhysics<3>;

} // namespace FeelModels
} // namespace Feel
